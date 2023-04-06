#ifndef _CGAL_IO_H
#define _CGAL_IO_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace fs = ::boost::filesystem;
// GEOMETRY //
typedef CGAL::Simple_cartesian<double>			Kernel;

// 3D //
typedef Kernel::Point_2               Point_2;
typedef Kernel::Point_3				        Point_3;
typedef Kernel::Vector_3				Vector_3;


// SURFACE_MESH //
typedef CGAL::Surface_mesh<Point_3>			Surface_mesh;
typedef Surface_mesh::Vertex_index				Vertex;
typedef Surface_mesh::Face_index				Face;
typedef Surface_mesh::Property_map<Vertex, Point_3>		VProp_geom;
typedef boost::graph_traits<Surface_mesh>::vertex_iterator vertex_iterator;

bool CGAL_read_OBJ(std::istream& is, Surface_mesh& sm);
bool outfileExists(fs::path outFilePath, const int size, std::string printDesc);
bool infileExists(fs::path inFilePath, const int size, std::string errDesc);
std::vector<std::string> splitString(std::string str, std::string delimiters);
bool meshLoader(fs::path meshFile, Surface_mesh& loadedMesh, std::string fileDesc, bool bdebug = true);
bool saveMesh(fs::path meshFile, Surface_mesh& sm, std::string fileDesc, std::stringstream& LogFile);

#endif