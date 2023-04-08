#ifndef _CGAL_CUT_UTILS_H
#define _CGAL_CUT_UTILS_H
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include "cgal_include.h"
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
//all the type we need for dynamic changing the cut path
typedef typename Kernel::FT                                                    NT;
typedef CGAL::dynamic_face_property_t<NT>                                      Face_NT_tag;
typedef typename boost::property_map<Seam_mesh, Face_NT_tag>::type             Face_NT_map;

typedef std::map<SM_vertex_descriptor, int>                              VertexIndexMap;
typedef boost::associative_property_map<VertexIndexMap>               VertexIdPropertyMap;
typedef boost::graph_traits<Surface_mesh>::vertex_iterator           SM_vertex_iterator;
namespace SMP = CGAL::Surface_mesh_parameterization;
namespace PMP = CGAL::Polygon_mesh_processing;
typedef typename SMP::internal::Kernel_traits<Seam_mesh>::PPM                       PPM;
typedef typename boost::property_traits<PPM>::reference                        PPM_ref;
typedef std::unordered_set<vertex_descriptor>                                  Vertex_set;
static Face_NT_map m_face_areas;

void IndexConvert2Descriptor(
	Surface_mesh&tmesh,
	std::vector<std::vector<int>>&init_cutpaths,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths);


void update_cutpaths(
	Seam_mesh&tmesh,
	Face_NT_map& face_L2_map,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths
);


Face_NT_map get_dist_map(
	Seam_mesh& seam_mesh,
	UV_pmap& uv_map
);

NT get_a(
	const std::array<const Point_2*, 3>& uv_points
);

Point_3 get_ss(
	const std::array<const Point_3*, 3>& mesh_points,
	const std::array<const Point_2*, 3>& uv_points,
	const NT den
);


Point_3 get_st(
	const std::array<const Point_3*, 3>& mesh_points,
	const std::array<const Point_2*, 3>& uv_points,
	const NT den
);

NT inner_product(
	const Point_3& pointA,
	const Point_3& pointB
);

double get_error(
	Seam_mesh& seam_mesh,
	UV_pmap& uv_map,
	halfedge_descriptor& bhd
);

NT initialize_faces_areas(
	const std::vector<face_descriptor>& face_range,
	Seam_mesh& tmesh
);


NT compute_area_distortion(
	const std::vector<face_descriptor>& face_range,
	const NT A_3D,
	Seam_mesh& tmesh,
	UV_pmap& uvmap
);

NT
compute_stretch(
	Seam_mesh& tmesh, 
	Face_NT_map& face_L2_map
);
#endif