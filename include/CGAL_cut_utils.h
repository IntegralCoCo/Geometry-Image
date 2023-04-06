#ifndef _CGAL_CUT_UTILS_H
#define _CGAL_CUT_UTILS_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include "cgal_include.h"
//all the type we need for dynamic changing the cut path
typedef typename Kernel::FT                                                    NT;
typedef CGAL::dynamic_face_property_t<NT>                                      Face_NT_tag;
typedef typename boost::property_map<Seam_mesh, Face_NT_tag>::type             Face_NT_map;

typedef std::map<SM_vertex_descriptor, int>                              VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
typedef boost::graph_traits<Surface_mesh>::vertex_iterator           SM_vertex_iterator;

void IndexConvert2Descriptor(
	Surface_mesh&tmesh,
	std::vector<std::vector<int>>&init_cutpaths,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths);


void update_cutpaths(
	Seam_mesh&tmesh,
	Face_NT_map& face_L2_map,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths
);




#endif