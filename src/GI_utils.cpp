#include "GI_utils.h"
#include <Eigen/Geometry>
#include <CGAL/license/Surface_mesh_parameterization.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/circulator.h>

#include <boost/iterator/function_output_iterator.hpp>
int Ceil(cv::Mat A) {
    int tmp = ceil(A.at<float>(0, 0));
    for (int i = 1; i < A.cols; ++i) {
        tmp = ceil(A.at<float>(0, i)) > tmp ? ceil(A.at<float>(0, i)) : tmp;
    }
    return tmp;
}

int Floor(cv::Mat A) {
    int tmp = floor(A.at<float>(0, 0));
    for (int i = 1; i < A.cols; ++i) {
        tmp = floor(A.at<float>(0, i)) < tmp ? floor(A.at<float>(0, i)) : tmp;
    }
    return tmp;
}

static Eigen::Vector3d convertToBarycentricCoordinates2D(std::vector<Point_2>&face_P ,double u,double v)
{
    Eigen::Vector3d result;
    Eigen::Vector3d edge1, edge2, p;
    edge1(0) = (face_P[1] - face_P[0])[0];
    edge1(1) = (face_P[1] - face_P[0])[1];
    edge1(2) = 0;

    edge2(0) = (face_P[2] - face_P[0])[0];
    edge2(1) = (face_P[2] - face_P[0])[1];
    edge2(2) = 0;

    p(0) = (u - face_P[0][0]);
    p(1) = (v - face_P[0][1]);
    p(2) = 0;

    double d00, d01, d11, d20, d21,denom;
    d00 = edge1.dot(edge1);
    d01 = edge1.dot(edge2);
    d11 = edge2.dot(edge2);
    d20 = p.dot(edge1);
    d21 = p.dot(edge2);
    denom = d00 * d11 - d01 * d01;

    result(1) = (d11 * d20 - d01 * d21) / denom;
    result(2) = (d00 * d21 - d01 * d20) / denom;
    result(0) = 1.0f - result(1) - result(2);

    return result;
}
bool PInsideTriangle(double u,double v, std::vector<Point_2>&face_P)
{
    const float EPSILON = 0.0000001f;
    Eigen::Vector3d edge1, edge2, edge3, h, s, q;

    edge1(0) = (face_P[1] - face_P[0])[0];
    edge1(1) = (face_P[1] - face_P[0])[1];
    edge1(2) = 0;

    edge2(0) = (face_P[2] - face_P[1])[0];
    edge2(1) = (face_P[2] - face_P[1])[1];
    edge2(2) = 0;

    edge3(0) = (face_P[0] - face_P[2])[0];
    edge3(1) = (face_P[0] - face_P[2])[1];
    edge3(2) = 0;

    h(0) = u - face_P[0][0];
    h(1) = v - face_P[0][1];
    h(2) = 0;

    s(0) = u - face_P[1][0];
    s(1) = v - face_P[1][1];
    s(2) = 0;

    q(0) = u - face_P[2][0];
    q(1) = v - face_P[2][1];
    q(2) = 0;


    double zh, zs, zq;
    zh = edge1.cross(h)[2];
    zs = edge2.cross(s)[2];
    zq = edge3.cross(q)[2];

    if ((zh >= EPSILON && zs >= EPSILON && zq >= EPSILON) || (zh <= -EPSILON && zs <= -EPSILON && zq <= -EPSILON))
    {
        return true;
    }
    else {
        return false;
    }
}

bool mesh2GI_planeSur(Seam_mesh& seam_mesh, halfedge_descriptor& bhd, UV_pmap& uvmap, const std::string OutGIPath, int im_size)
{
    boost::vertex_point_t VPT = boost::vertex_point;
    auto seam_v_map = CGAL::get(VPT, seam_mesh);
    typedef std::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);

    std::unordered_set<vertex_descriptor> vertices;
    std::vector<face_descriptor> faces;

    CGAL::Surface_mesh_parameterization::internal::Containers_filler<Seam_mesh> fc(seam_mesh, vertices, &faces);
    CGAL::Polygon_mesh_processing::connected_component(
        face(opposite(bhd, seam_mesh), seam_mesh),
        seam_mesh,
        boost::make_function_output_iterator(fc));


    // matrices for GI
    cv::Mat M0 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    cv::Mat M1 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    cv::Mat M2 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    std::map<int, cv::Mat> outputMap;
    outputMap[0] = M0;
    outputMap[1] = M1;
    outputMap[2] = M2;

    double u, v;
    for (int y = 0; y < im_size; y++)
    {
        v = (double) y / (im_size - 1);
        printf("rendering:%d/%d\n", y, im_size);
        for (int x = 0; x < im_size; x++)
        {
            u = (double) x / (im_size - 1);
            //judge which triangle this pixel is located at
            for (face_descriptor fd : faces) {
                std::vector<Point_2> face_v;
                std::vector<Point_3> face_v3D;
                halfedge_descriptor hd = halfedge(fd, seam_mesh);
                int vt_count = 0;
                for (vertex_descriptor vd : vertices_around_face(hd, seam_mesh)) {
                    Point_2 P2D = get(uvmap, vd);
                    Point_3 P3D = get(seam_v_map, vd);
                    face_v.push_back(P2D);
                    face_v3D.push_back(P3D);
                }
                if (face_v.size() < 3)
                {
                    break;
                }
                if (PInsideTriangle(u, v, face_v))
                {
                    Eigen::Vector3d baycoord = convertToBarycentricCoordinates2D(face_v, u,v);
                    //x
                    outputMap[0].at<float>(y, x) = face_v3D[0][0] * baycoord(0) +
                                                   face_v3D[1][0] * baycoord(1)+
                                                   face_v3D[2][0] * baycoord(2);
                    outputMap[1].at<float>(y, x) = face_v3D[0][1] * baycoord(0) +
                                                   face_v3D[1][1] * baycoord(1) +
                                                   face_v3D[2][1] * baycoord(2);
                    outputMap[2].at<float>(y, x) = face_v3D[0][2] * baycoord(0) +
                                                   face_v3D[1][2] * baycoord(1) +
                                                   face_v3D[2][2] * baycoord(2);
                    break;
                };
            }
        }
    }

    combineNSave(outputMap, OutGIPath, ",save to GI", im_size);
    return true;
}


void filterOutputMap(cv::Mat& A, cv::Mat& mask_value, cv::Mat& mask_NaN,
    cv::Mat& kernel, int nIter) {
    // compute mean of value outputMap and fed those values to NaN
    float mean_outValue = cv::mean(A, mask_value)[0];
    cv::Mat out_tmp(A.size(), A.type(), cv::Scalar(mean_outValue));
    out_tmp.copyTo(A, mask_NaN);
    out_tmp.release();

    // Perform convolution using a ones filter
    cv::Mat outputMapTmp;
    A.copyTo(outputMapTmp);
    for (int i = 0; i < nIter; i++) {
        cv::filter2D(outputMapTmp, outputMapTmp, -1, kernel, cv::Point(-1, -1),
            0, cv::BORDER_DEFAULT);
        A.copyTo(outputMapTmp, mask_value);
    }
    outputMapTmp.copyTo(A);
    outputMapTmp.release();
}


void combineNSave(std::map<int, cv::Mat>& outMap, std::string meshFileFlatGI, std::string desc,int im_size) {
    // this check is to ensure any of the previous files are not overwritten
    //if (outfileExists(meshFileFlatGI, 10, desc + "ED"))
    //    return;

    // statistics for the three channels
    double minVal[3];
    double maxVal[3];
    // for each Dimension of 3D mesh
    for (int dim = 0; dim < 3; ++dim) {
        // find the min and max in each dimension/plane seperately
        cv::minMaxLoc(outMap[dim], &minVal[dim], &maxVal[dim], NULL, NULL);
        // subtract with the min so that the new min is Zero, equivalent to translation in 3D
        cv::subtract(outMap[dim], minVal[dim], outMap[dim]);
    }

    cv::Mat in[] = { outMap[2], outMap[1], outMap[0] };
    int from_to[] = { 0, 0, 1, 1, 2, 2 };
    cv::Mat M = cv::Mat::zeros(im_size, im_size, CV_32FC3); // 3D Geometry Image
    cv::mixChannels(in, 3, &M, 1, from_to, 3);

    // divide by maximum so that new maximum is one
    // equivalent to uniform scaling as it is applied to all the dimensions
    cv::divide(M, newMax(minVal, maxVal), M);

    cv::Mat MM(im_size, im_size, CV_8UC3);
    M.convertTo(MM, CV_8UC3, 255);
    // scale to 255 is required or else image wont be visible in other viewer

    // swap color channels because opencv saves as BGR
    cvtColor(MM, MM, cv::COLOR_RGB2BGR);
    //std::vector<int> compression_params;
    //compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    //compression_params.push_back(0);

    cv::imwrite(meshFileFlatGI, MM);
    std::cout << desc << std::flush;
}

double newMax(double minVal[3], double maxVal[3]) {
    double tmp, tmpVal[3];
    int J;
    for (int j = 0; j < 3; j++) {
        tmpVal[j] = maxVal[j] - minVal[j];
        if (true && j == 0)
            // if its a slice means that the x axis should be giving double length
            tmpVal[j] *= 2;
        if (j == 0) {
            tmp = tmpVal[j];
            J = 0;
        }
        else if (tmpVal[j] > tmp) {
            tmp = tmpVal[j];
            J = j;
        }
    }
    return tmpVal[J];
}

bool realmesh2GI(Seam_mesh& seam_mesh, halfedge_descriptor& bhd, UV_pmap& uvmap, const std::string OutGIPath, int im_size)
{
    boost::vertex_point_t VPT = boost::vertex_point;
    auto seam_v_map = CGAL::get(VPT, seam_mesh);
    typedef std::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);

    std::unordered_set<vertex_descriptor> vertices;
    std::vector<face_descriptor> faces;

    CGAL::Surface_mesh_parameterization::internal::Containers_filler<Seam_mesh> fc(seam_mesh, vertices, &faces);
    CGAL::Polygon_mesh_processing::connected_component(
        face(opposite(bhd, seam_mesh), seam_mesh),
        seam_mesh,
        boost::make_function_output_iterator(fc));

    //std::ostringstream out_vertices, out_faces;
    //std::size_t vertices_counter = 0, faces_counter = 0;
    cv::Mat M0 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    cv::Mat M1 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    cv::Mat M2 = cv::Mat::zeros(im_size, im_size, CV_32FC1);  // Geometry Image
    std::map<int, cv::Mat> outputMap;
    outputMap[0] = M0;
    outputMap[1] = M1;
    outputMap[2] = M2;

    cv::Mat P = cv::Mat::zeros(2, 3, CV_32F); // Current 2D point
    cv::Mat V = cv::Mat::zeros(3, 3, CV_32F); // current dimension value of 3D Mesh
    cv::Mat nV = cv::Mat::zeros(3, 3, CV_32F);  // current dimension value of 3D Mesh Normals
    int f_idx = 0;
    cv::Mat Mnb = cv::Mat::zeros(im_size, im_size, CV_32FC1); // Geometry Image # pts calculator
    for (face_descriptor fd : faces) {
        halfedge_descriptor hd = halfedge(fd, seam_mesh);
        int vt_count = 0;
        for (vertex_descriptor vd : vertices_around_face(hd, seam_mesh)) {
            Point_2 P2D = get(uvmap, vd);
            Point_3 P3D = get(seam_v_map, vd);
            P.at<float>(0, vt_count) = P2D[0] * (im_size - 1);
            P.at<float>(1, vt_count) = P2D[1] * (im_size - 1);
            for (int dim = 0; dim < 3; ++dim) {
                //each Dimension of 3D mesh
                V.at<float>(dim, vt_count) = P3D[dim];
                //if (useNormal)
                //    nV.at<float>(dim, vt_count) = Mesh_3D_nm[vd_3D][dim];
            }
            vt_count++;
        }
        cv::Mat pos, z1, z2;
        // Dimensions of the mesh grid
        int n_rows = Ceil(P.row(0)) - Floor(P.row(0)) + 1;
        int n_cols = Ceil(P.row(1)) - Floor(P.row(1)) + 1;
        if (n_rows <= 0 || n_cols <= 0)
            continue;
        try {
            z1.create(n_rows, n_cols, CV_32FC1);
            z2.create(n_rows, n_cols, CV_32FC1);
        }
        catch (...) {
            std::cerr << "Couldn't create z1 or z2\n";
            return false;
        }
        for (int i = 0; i < z1.rows; ++i) {
            for (int j = 0; j < z1.cols; ++j) {
                z1.at<float>(i, j) = Floor(P.row(0)) + i;
                z2.at<float>(i, j) = Floor(P.row(1)) + j;
            }
        }

        z1 = z1.t();
        z1 = z1.reshape(0, 1);
        z2 = z2.t();
        z2 = z2.reshape(0, 1);
        cv::vconcat(z1, z2, pos);

        // barycentric coords computation
        cv::Mat tmp1, tmp2, c;
        cv::vconcat(cv::Mat::ones(1, P.cols, CV_32F), P, tmp1);
        cv::vconcat(cv::Mat::ones(1, pos.cols, CV_32F), pos, tmp2);
        cv::solve(tmp1, tmp2, c);

        // restrict the BC to inside triangle
        cv::Mat pos1, c1;
        c1.convertTo(c1, CV_32F);
        pos1.convertTo(pos1, CV_32F);
        for (int i = 0; i < c.cols; i++) {
            // cutoff was previously taken from octave but should be different for c++ beacuse of difference in datatypes
            if (c.at<float>(0, i) >= -0.000022204
                && c.at<float>(1, i) >= -0.000022204
                && c.at<float>(2, i) >= -0.000022204
                && !((c.at<float>(0, i) == 0) && (c.at<float>(1, i) == 0)
                    && (c.at<float>(2, i) == 0))) {
                cv::Mat trash1 = c.col(i).t();
                c1.push_back(trash1);
                cv::Mat trash2 = pos.col(i).t();
                pos1.push_back(trash2);
            }
        }
        if (c1.empty() || pos1.empty()) {
            pos1.release();
            c1.release();
            f_idx++;
            continue;
        }
        c = c1.t();
        pos = pos1.t();
        pos1.release();
        c1.release();

        // restrict the pos to inside the image
        // TODO: this condition never seems to happen
        for (int i = 0; i < pos.cols; i++) {
            if (pos.at<float>(0, i) >= 0) {
                if (pos.at<float>(0, i) < im_size) {
                    if (pos.at<float>(1, i) >= 0) {
                        if (pos.at<float>(1, i) < im_size) {
                            cv::Mat trash1 = c.col(i).t();
                            c1.push_back(trash1);
                            cv::Mat trash2 = pos.col(i).t();
                            pos1.push_back(trash2);

                        }
                    }
                }
            }
        }

        if (c1.empty() || pos1.empty()) {
            pos1.release();
            c1.release();
            f_idx++;
            continue;
        }

        c = c1.t();
        pos = pos1.t();
        pos1.release();
        c1.release();

        // Now the value assignment has to be done
        for (int i = 0; i < pos.cols; i++) {
            int r_idx = (int)pos.at<float>(0, i);
            int c_idx = (int)pos.at<float>(1, i);
            for (int dim = 0; dim < 3; ++dim) {  //each Dimension of 3D mesh
                outputMap[dim].at<float>(r_idx, c_idx) =
                    outputMap[dim].at<float>(r_idx, c_idx)
                    + V.at<float>(dim, 0) * c.at<float>(0, i)
                    + V.at<float>(dim, 1) * c.at<float>(1, i)
                    + V.at<float>(dim, 2) * c.at<float>(2, i);
            }
            Mnb.at<float>(r_idx, c_idx) = Mnb.at<float>(r_idx, c_idx) + 1;
        }
        f_idx++;
    }

    cv::Mat tmp_mask_value(im_size, im_size, CV_8UC1);
    cv::Mat tmp_mask_NaN(im_size, im_size, CV_8UC1);
    cv::threshold(Mnb, tmp_mask_value, 0, 255, cv::THRESH_BINARY);
    cv::threshold(Mnb, tmp_mask_NaN, 0, 255, cv::THRESH_BINARY_INV);

    cv::Mat mask_NaN, mask_value;
    tmp_mask_NaN.convertTo(mask_NaN, CV_8UC1);
    tmp_mask_value.convertTo(mask_value, CV_8UC1);
    tmp_mask_NaN.release();
    tmp_mask_value.release();

    // Specify convolution filter size required to compensate for no values in geometry image
    int filterX = 3;
    int filterY = 3;
    int nIter = 20;
    cv::Mat kernel = cv::Mat::ones(filterX, filterY, CV_32F) / (float)(filterX * filterY);

    for (int dim = 0; dim < 3; ++dim) {  //each Dimension of 3D mesh
        if (outputMap[dim].cols == 0) {
            std::cerr << "  OutputMap is empty" << std::endl;
            return false;
        }

        // divide by the number of vertices
        outputMap[dim] /= Mnb;
        // filter the NaNs
        // 1. for GI
        filterOutputMap(outputMap[dim], mask_value, mask_NaN, kernel, nIter);
        // 2. if required then for mormal GI
    } //each dimension

    // 1. for GI
    combineNSave(outputMap, OutGIPath, ", savedGI", im_size);
    // 2. if required for Normal GI
    return true;
}