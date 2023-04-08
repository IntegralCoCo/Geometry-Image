#include "square_border_parameterization.h"
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include "GI_utils.h"
#include "CGAL_cut_utils.h"
Parameterization::Parameterization(fs::path inputPath, fs::path OutGIPath, bool& useNormal, int& im) :
    useNormal(useNormal), im_size(im) {
    this->inputPath = inputPath;
    this->OutGIPath = (OutGIPath / inputPath.stem()).string() + "_flatGI.png";
    this->OutnGIPath = (OutGIPath.parent_path() / OutGIPath.stem()).string() + "_" + std::to_string(im_size) + "_nflatGI.png";


}

Parameterization::~Parameterization() {
    // TODO Auto-generated destructor stub
}

bool Parameterization::surfaceParameterizeIterative(std::vector<std::vector<int>>&init_cutpaths) {
    // check if the parameterised file already exist;
    // read input
    Surface_mesh sm;
    if (!meshLoader(inputPath, sm, " input mesh for parameterization"))
        return false;

    //define seam mesh
    Seam_edge_uhm seam_edge_uhm(false);
    Seam_edge_pmap seam_edge_pm(seam_edge_uhm);

    Seam_vertex_uhm seam_vertex_uhm(false);
    Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);
    Seam_mesh seam_mesh(sm, seam_edge_pm, seam_vertex_pm);
    

    std::vector<std::vector<SM_vertex_descriptor>> cutpaths;
    IndexConvert2Descriptor(sm, init_cutpaths, cutpaths);
    //boost::vertex_point_t VPT = boost::vertex_point;
    //auto sm_v_map = CGAL::get(VPT, sm);
    //for (SM_vertex_descriptor vd : sm.vertices())
    //{

    //    std::cout << vd.idx() << std::endl;
    //    break;
    //}
    //vertex_descriptor vd;
    //boost::vertex_point_t t;
    //auto p = CGAL::get(t,seam_mesh);
    //Point_3 tmp = get(p, vd);

    //add seam
    //auto first_path = cutpaths.front();
    //for (int i = 0; i < first_path.size() - 1; i++)
    //{
    //    if (!seam_mesh.add_seam(first_path[i], first_path[i + 1])) {
    //        std::cout << "WARNING: Deprecated add seam! Ignore it(This may cause serious error)" << std::endl;
    //        std::cout << first_path[i].id() << std::endl;
    //        std::cout << first_path[i + 1].id() << std::endl;
    //    };
    //}
    double last_stretch = std::numeric_limits<double>::max();
    halfedge_descriptor output_bhd;
    UV_uhm  output_uv_uhm;
    UV_pmap output_uv_map;
    while (1)
    {
        //SM_halfedge_descriptor smhd;
        auto last_path = cutpaths.back();
        for (int i = 0; i < last_path.size() - 1; i++)
        {
            if (!seam_mesh.add_seam(last_path[i], last_path[i + 1])) {
                std::cout << "WARNING: Deprecated add seam! Ignore it(This may cause serious error)" << std::endl;
                std::cout << last_path[i].idx() << std::endl;
                std::cout << last_path[i+1].idx() << std::endl;
            };
        }
        //smhd = seam_mesh.add_seams(last_path.begin(), last_path.end());
        //if (smhd == SM_halfedge_descriptor()) {
        //    std::cerr << "Warning: No seams in input" << std::endl;
        //}
        halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(seam_mesh).first;
        UV_uhm  uv_uhm;
        UV_pmap uv_map(uv_uhm);
        SMP::Error_code err;
        std::cout << "Circular Parameterization started" << std::endl;
        Circular_Parameterizer c_param;
        try {
            err = SMP::parameterize(seam_mesh, c_param, bhd, uv_map);
        }
        catch (...) {
            std::cerr << "  SMP::parameterize didn't succeed\n";
            return false;
        }
        std::cout << "Circular Parameterization finished" << std::endl;
        if (err != SMP::OK) {
            std::cerr << "  Error: " << SMP::get_error_message(err) << "\n";
            return false;
        }
        Face_NT_map face_dist_map = get_dist_map(seam_mesh, uv_map);
        update_cutpaths(seam_mesh, face_dist_map, cutpaths);
        
        halfedge_descriptor temp_bhd = CGAL::Polygon_mesh_processing::longest_border(seam_mesh).first;
        // final square border parameterization
        // The 2D points of the uv parameterization will be written into this map
        UV_uhm  temp_uv_uhm;
        UV_pmap temp_uv_map(temp_uv_uhm);

        std::cout << "Square Parameterization started" << std::endl;
        try {
            err = Authalic_Parameterizer().parameterize(seam_mesh, temp_bhd, temp_uv_map);
        }
        catch (...) {
            std::cerr << "  SMP::parameterize didn't succeed\n";
            return false;
        }
        std::cout << "Parameterization finished" << std::endl;
        if (err != SMP::OK) {
            std::cerr << "  Error: " << SMP::get_error_message(err) << "\n";
            return false;
        }
        Face_NT_map s_face_dist_map = get_dist_map(seam_mesh, temp_uv_map);
        double cur_stretch = get_error(seam_mesh, temp_uv_map, temp_bhd);
        //double cur_stretch = compute_stretch(seam_mesh, s_face_dist_map);
        if (last_stretch < cur_stretch)
        {
            output_uv_map = UV_pmap(output_uv_uhm);
            std::cout << "stretch rises, end the" << std::endl;
            break;
        }
        output_uv_uhm = temp_uv_uhm;
        output_bhd = temp_bhd;
    	std::cout << "Last stretch:" << last_stretch << " " << "Current stretch:" << cur_stretch << std::endl;
    	last_stretch = cur_stretch;
    }

    std::ofstream out("result.off");
    SMP::IO::output_uvmap_to_off(seam_mesh, output_bhd, output_uv_map, out);
    realmesh2GI(seam_mesh, output_bhd, output_uv_map, "GI_test.bmp", 257);

    // reconstruct mesh from uv map


    return true;
}