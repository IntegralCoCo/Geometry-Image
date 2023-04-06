#ifndef _CGAL_INCLUDE_H
#define _CGAL_INCLUDE_H

#include <array>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <memory>
#include <random>
#include <signal.h>
#include <stdexcept>
#include <stdlib.h>
#include <string>
#include <sys/stat.h>
#include <vector>


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include "cgal_IO.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Plane_3 K_Plane_3;

typedef boost::graph_traits<Surface_mesh>::vertex_descriptor SM_vertex_descriptor;
typedef boost::graph_traits<Surface_mesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<Surface_mesh>::face_descriptor SM_face_descriptor;
typedef boost::graph_traits<Surface_mesh>::halfedge_descriptor SM_halfedge_descriptor;

//typedef Surface_mesh::Property_map<SM_vertex_descriptor, int> SM_vimap;



typedef CGAL::Unique_hash_map<SM_halfedge_descriptor, Point_2> UV_uhm;
typedef boost::associative_property_map<UV_uhm> UV_pmap;

typedef CGAL::Unique_hash_map<SM_edge_descriptor, bool> Seam_edge_uhm;
typedef CGAL::Unique_hash_map<SM_vertex_descriptor, bool> Seam_vertex_uhm;
typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm>Seam_vertex_pmap;
//
//
typedef CGAL::Seam_mesh<Surface_mesh, Seam_edge_pmap, Seam_vertex_pmap>     Seam_mesh;
//
typedef boost::graph_traits<Seam_mesh>::vertex_descriptor         vertex_descriptor;
typedef boost::graph_traits<Seam_mesh>::halfedge_descriptor       halfedge_descriptor;
typedef boost::graph_traits<Seam_mesh>::edge_descriptor           edge_descriptor;
typedef boost::graph_traits<Seam_mesh>::face_descriptor           face_descriptor;


#endif