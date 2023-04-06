#include <igl/harmonic.h>
#include <igl/MappingEnergyType.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/slim.h>
#include <igl/arap.h>

void param_2d_SLIM(const Eigen::MatrixXd&V, 
				   const Eigen::MatrixXi&F, 
				   const Eigen::MatrixXd&uv_init,
				   igl::SLIMData& sData,
	               int max_count,
				   bool isDebug = false) {

	//parametrization
	sData.slim_energy = igl::MappingEnergyType::SYMMETRIC_DIRICHLET;
	Eigen::VectorXi b;
	Eigen::VectorXd bc;
	igl::slim_precompute(V, F, uv_init, sData, igl::MappingEnergyType::SYMMETRIC_DIRICHLET,b,bc,0);
	igl::slim_solve(sData, max_count);
	if (isDebug)
	{
		printf("it times:%d, current energy:%.5f\n", max_count, sData.energy / sData.mesh_area);
	}
}


//void param_2d_ARAP(const Eigen::MatrixXd& V,
//				   const Eigen::MatrixXi& F,
//				   const Eigen::MatrixXd& uv_init,
//				   igl::ARAPData& arap_data,
//				   int max_count,
//	               bool isDebug = false)
//{
//	arap_data.with_dynamics = true;
//	Eigen::VectorXi b = Eigen::VectorXi::Zero(0);
//	Eigen::MatrixXd bc = Eigen::MatrixXd::Zero(0, 0);
//
//	arap_data.max_iter = max_count;
//	// 2 means that we're going to *solve* in 2d
//	arap_precomputation(V, F, 2, b, arap_data);
//	igl::arap_solve(bc, arap_data, uv_init);
//	if (isDebug)
//	{
//		printf("it times:%d, current energy:%.5f\n", max_count, arap_data.energy);
//	}
//}


void iterativelyUpdateCutPaths(std::string mesh_path, 
	                           std::vector<std::vector<int>>&init_cutpaths,
							   std::string parameterizer
)
{
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;
	Eigen::MatrixXd Vn;
	Eigen::MatrixXi Fn;
	Eigen::MatrixXd init_uvmap;
	bool isgenus_zero = false;
	if (igl::readOFF(mesh_path, V, F))
	{
		std::cout << "loaded mesh" << std::endl;
		printf("mesh Vertexs:%ld\n", V.rows());
		printf("mesh Faces:%ld\n", F.rows());
	}
	//cutpath
	igl::cut_to_disk(F, init_cutpaths);
	//if init_cutpath returns nothing,then we randomly select a vertex and it's two adjacent edge(or 1?)
	srand(time(0));
	if (init_cutpaths.size() == 0)
	{
		int rand_Face = rand() % F.size() + 1;
		std::vector<int> rand_path{ F(rand_Face,0),F(rand_Face,1),F(rand_Face,2) };
		init_cutpaths.push_back(rand_path);
		isgenus_zero = true;
	}
	const size_t num_faces = F.rows();
	Eigen::MatrixXi cut_mask(num_faces, 3);
	cut_mask.setZero();
	cut_utils::assert_is_disk(V, F, init_cutpaths, cut_mask);
	igl::cut_mesh(V, F, cut_mask, Vn, Fn);
	//parametrization
	Eigen::VectorXi bnd;
	igl::boundary_loop(Fn, bnd);
	Eigen::MatrixXd bnd_uv;
	igl::map_vertices_to_circle(Vn, bnd, bnd_uv);
	igl::harmonic(Vn, Fn, bnd, bnd_uv, 1, init_uvmap);
	igl::SLIMData sData;
	igl::ARAPData aData;
	int max_count = 5;
	if (parameterizer == "slim")
	{
		param_2d_SLIM(Vn, Fn, init_uvmap, sData, max_count, true);
		std::cout << "Finish parameterization" << std::endl;
		std::cout << "current energy:" << sData.energy << std::endl;
	}
	else if (parameterizer == "arap")
	{
		//param_2d_ARAP(Vn, Fn, init_uvmap, aData, max_count, true);
		std::cout << "Finish parameterization" << std::endl;
		std::cout << "current energy:" << aData.energy << std::endl;
	}

	//we now update cut path base on distortion
	int update_maxcount = 1;
	std::cout << "iterative update cut path" << std::endl;
	int cur = 0;

	while (init_cutpaths[init_cutpaths.size() - 1].size() > 3 || isgenus_zero)
	{
		isgenus_zero = false;
		cut_utils::updata_cut_path(V, F, init_cutpaths, sData);
		cut_mask.setZero();
		cut_utils::assert_is_disk(V, F, init_cutpaths, cut_mask);
		igl::cut_mesh(V, F, cut_mask, Vn, Fn);

		igl::SLIMData sData_update;
		igl::ARAPData aData_update;

		igl::boundary_loop(Fn, bnd);
		igl::map_vertices_to_circle(Vn, bnd, bnd_uv);
		igl::harmonic(Vn, Fn, bnd, bnd_uv, 1, init_uvmap);

		if (parameterizer == "slim")
		{
			param_2d_SLIM(Vn, Fn, init_uvmap, sData_update, max_count, true);
			std::cout << "Finish parameterization" << std::endl;
			std::cout << "current energy:" << sData_update.energy << std::endl;
		}
		else if (parameterizer == "arap")
		{
			//param_2d_ARAP(Vn, Fn, init_uvmap, aData, max_count, true);
			std::cout << "Finish parameterization" << std::endl;
			std::cout << "current energy:" << aData_update.energy << std::endl;
		}

		//show out
		cur++;
		if (cur > update_maxcount)
		{
			break;
		}
	}
}