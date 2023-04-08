#ifndef _SQUARE_BORDER_H
#define _SQUARE_BORDER_H

#include "cgal_include.h"
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/surface_mesh_parameterization/Circular_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Mean_value_coordinates_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Iterative_authalic_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>


#include<vector>



namespace SMP = CGAL::Surface_mesh_parameterization;
typedef SMP::Square_border_arc_length_parameterizer_3<Seam_mesh> Square_Border_parameterizer;
typedef SMP::Mean_value_coordinates_parameterizer_3<Seam_mesh,Square_Border_parameterizer> Square_Parameterizer;
typedef SMP::Iterative_authalic_parameterizer_3<Seam_mesh> Authalic_Parameterizer;
typedef SMP::Circular_border_arc_length_parameterizer_3<Seam_mesh> Circular_Border_parameterizer;
typedef SMP::Mean_value_coordinates_parameterizer_3<Seam_mesh, Circular_Border_parameterizer>Circular_Parameterizer;


namespace PMP = CGAL::Polygon_mesh_processing;



class Parameterization {
public:
	Parameterization(fs::path inputPath, fs::path outGIPath, bool& useNormal, int& im);
	virtual ~Parameterization();
	bool surfaceParameterizeIterative(std::vector<std::vector<int>>& init_cutpaths);
	bool mesh2GI();
	bool GI2off();


private:
	double newMax(double minVal[3], double maxVal[3]);
	bool readGI(int downScaleFactor = 1);
	bool addVerticestoSM(Surface_mesh& sm);

	fs::path inputPath; // input path
	bool useNormal; // use normals with geometry image
	int im_size;

	fs::path GIPath; // surface paramterized mesh
	std::string OutGIPath;
	std::string OutnGIPath;  // normal encoded geometry image

	std::stringstream verticesSS, normalSS, faceSS; //strings to read mesh from GI
};


#endif