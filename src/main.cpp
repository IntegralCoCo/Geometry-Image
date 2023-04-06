#include<igl/cut_to_disk.h>
#include<igl/cut_mesh.h>
#include<igl/readOFF.h>
#include<Eigen/Dense>
#include<fstream>
#include<string>
#include<vector>
#include<assert.h>
#include"square_border_parameterization.h"
#include"GI_utils.h"
#include"CGAL_cut_utils.h"
#include<stdlib.h>
#include<opencv2/opencv.hpp>

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd Vn;
Eigen::MatrixXi Fn;
Eigen::MatrixXd init_uvmap;
//void writeUVmap(const std::string& path, Eigen::MatrixXd& uvmap, Eigen::MatrixXi& F);
//void writeSelectiontxt(const std::vector<std::vector<int>>& cutpaths);//the file path must be settled down
int main(int argc, char* argv[]) {
	bool isgenus_zero = false;
	std::string mesh_path = "bunny.off";
	if (igl::readOFF(mesh_path, V, F))
	{
		std::cout << "loaded mesh" << std::endl;
		printf("mesh Vertexs:%ld\n", V.rows());
		printf("mesh Faces:%ld\n", F.rows());
	}
	//cutpath
	std::vector<std::vector<int>> init_cutpaths;
	igl::cut_to_disk(F, init_cutpaths);
	srand(time(0));
	if (init_cutpaths.size() == 0)
	{
		int rand_fd = rand() % F.rows();
		init_cutpaths.push_back(std::vector<int>{F(rand_fd, 0), F(rand_fd, 1), F(rand_fd, 2)});
		//init_cutpaths.push_back(std::vector<int>{F(rand_fd, 1), F(rand_fd, 2)});
	}
	////if init_cutpath returns nothing,then we randomly select a vertex and it's two adjacent edge(or 1?)
	//srand(time(0));
	//if (init_cutpaths.size() == 0)
	//{
	//	int rand_Face = rand() % F.size() + 1;
	//	std::vector<int> rand_path{ F(rand_Face,0),F(rand_Face,1),F(rand_Face,2)};
	//	init_cutpaths.push_back(rand_path);
	//	isgenus_zero = true;
	//}
	//const size_t num_faces = F.rows();
	//Eigen::MatrixXi cut_mask(num_faces,3);
	//cut_mask.setZero();
	//cut_utils::assert_is_disk(V, F, init_cutpaths, cut_mask);
	//igl::cut_mesh(V, F, cut_mask, Vn, Fn);
	////parametrization
	//Eigen::VectorXi bnd;
	//igl::boundary_loop(Fn, bnd);
	//Eigen::MatrixXd bnd_uv;
	//igl::map_vertices_to_circle(Vn, bnd, bnd_uv);
	//igl::harmonic(Vn, Fn, bnd, bnd_uv, 1, init_uvmap);
	//igl::SLIMData sData;
	//int max_count = 1;
	//param_2d_SLIM(Vn, Fn,init_uvmap,sData, max_count,true);
	//writeUVmap("uvmap_init.obj", sData.V_o, Fn);
	////show out
	//std::cout << "Finish parameterization" << std::endl;
	//std::cout << "current energy:" << sData.energy << std::endl;
	////we now update cut path base on distortion
	//int update_maxcount = 1;
	//std::cout << "iterative update cut path" << std::endl;
	//int cur = 0;
	//while (init_cutpaths[init_cutpaths.size()-1].size() > 3 || isgenus_zero)
	//{
	//	isgenus_zero = false;
	//	cut_utils::updata_cut_path(V, F, init_cutpaths, sData);
	//	cut_mask.setZero();
	//	cut_utils::assert_is_disk(V, F, init_cutpaths, cut_mask);
	//	igl::cut_mesh(V, F, cut_mask, Vn, Fn);
	//	igl::SLIMData sData_update;
	//	igl::boundary_loop(Fn, bnd);
	//	igl::map_vertices_to_circle(Vn, bnd, bnd_uv);
	//	igl::harmonic(Vn, Fn, bnd, bnd_uv, 1, init_uvmap);
	//	param_2d_SLIM(Vn, Fn, init_uvmap, sData_update, max_count,true);
	//	writeUVmap("uvmap_update.obj", sData_update.V_o, Fn);
	//	//show out
	//	std::cout << "Finish parameterization" << std::endl;
	//	std::cout << "current energy:" << sData_update.energy << std::endl;
	//	cur++;
	//	if (cur > update_maxcount)
	//	{
	//		break;
	//	}
	//}
	//writeSelectiontxt(init_cutpaths);
	std::string outModelPath = "./";
	fs::path ModelFilePath("bunny.off");
	fs::path ParamFilePath("result.off");
	std::string outGIPath = "GI_test.png";
	bool normal = false;
	int im_size = 63;
	Parameterization PM(ModelFilePath, outGIPath,normal,im_size);
	if (!PM.surfaceParameterizeIterative(init_cutpaths))
	{
		return -1;
	};
	return 0;
}

//void writeUVmap(const std::string& path, Eigen::MatrixXd& uvmap, Eigen::MatrixXi& F)
//{
//	Eigen::MatrixXd V_uvmap(Vn.rows(), 3);
//	V_uvmap.setZero();
//	for (int i = 0; i < uvmap.rows(); i++)
//	{
//		for (int j = 0; j < 2; j++)
//		{
//			V_uvmap(i, j) = uvmap(i, j);
//		}
//	}
//	igl::writeOBJ(path, V_uvmap, Fn);
//}
//
//void writeSelectiontxt(const std::vector<std::vector<int>>& cutpaths)
//{
//	const std::string filepath = "cutpath.selection.txt";
//	std::fstream wdes(filepath, std::ios::out);
//	
//	for (auto path : cutpaths)
//	{
//		char pair[20] = {'\0'};
//		for (int i = 0; i < path.size() - 1; i++)
//		{
//			sprintf(pair, "%d %d ", path[i], path[i + 1]);
//			wdes << pair;
//		}
//		wdes << std::endl;
//	}
//}

