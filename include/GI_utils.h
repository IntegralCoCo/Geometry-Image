#ifndef _GI_UTILS_H
#define _GI_UTILS_H
#include <opencv2/opencv.hpp>
#include "cgal_IO.h"
#include "cgal_include.h"
int Ceil(cv::Mat A);
int Floor(cv::Mat A);
bool mesh2GI(fs::path meshPath, fs::path uvmapPath, const std::string OutGIPath, int im_size=63);
void filterOutputMap(cv::Mat& A, cv::Mat& mask_value, cv::Mat& mask_NaN,
    cv::Mat& kernel, int nIter);
double newMax(double minVal[3], double maxVal[3]);
void combineNSave(std::map<int, cv::Mat>& outMap, std::string meshFileFlatGI, std::string desc, int im_size); bool realmesh2GI(Seam_mesh& seam_mesh, halfedge_descriptor& bhd, UV_pmap& uvmap, const std::string OutGIPath, int im_size);
bool mesh2GI_planeSur(Seam_mesh& seam_mesh, halfedge_descriptor& bhd, UV_pmap& uvmap, const std::string OutGIPath, int im_size);
static Eigen::Vector3d convertToBarycentricCoordinates2D(std::vector<Point_2>& face_P, double u, double v);
bool PInsideTriangle(double u, double v, std::vector<Point_2>& face_P);
#endif