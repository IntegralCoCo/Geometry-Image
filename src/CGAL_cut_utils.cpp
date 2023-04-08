#include "CGAL_cut_utils.h"


bool isInsideCutpaths(
	const Surface_mesh&sm,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths,
	face_descriptor fd)
{
	SM_face_descriptor sm_fd = SM_face_descriptor(fd.idx());
	//randomly select a vertex in face_maxdist
	SM_halfedge_descriptor sm_hd = halfedge(sm_fd, sm);
	// the source point can not be inside the cutpaths
	// means all 
	for (auto temp : CGAL::vertices_around_face(sm_hd, sm))
	{
		
		for (auto line : cutpaths)
		{
			for (auto vd : line)
			{
				if (temp == vd)
				{
					return true;
				}
			}
		}
	}
	return false;
}
void IndexConvert2Descriptor(
	Surface_mesh& tmesh,
	std::vector<std::vector<int>>& init_cutpaths,
	std::vector<std::vector<SM_vertex_descriptor>>&cutpaths)
{
	
	typedef std::unordered_map<std::size_t, SM_vertex_descriptor> Index_Vertex_map;
	Index_Vertex_map ivum;
	boost::associative_property_map<Index_Vertex_map> ivmap(ivum);
	std::size_t vertices_counter = 0;
	for (SM_vertex_descriptor vd : tmesh.vertices())
	{
		put(ivmap,vertices_counter++, vd);
	}
	for (auto row : init_cutpaths)
	{
		std::vector<SM_vertex_descriptor> tmp;
		for (int idx : row)
		{
			tmp.push_back(get(ivmap, idx));
		}
		cutpaths.emplace_back(tmp);
	}
	
}

void update_cutpaths(
	Seam_mesh& tmesh,
	Face_NT_map& face_L2_map,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths
)
{
	const Surface_mesh& sm = tmesh.mesh();
	face_descriptor max_distFd = *tmesh.faces_begin();
	NT max_dist = get(face_L2_map, max_distFd);

	for (face_descriptor fd : tmesh.faces())
	{
		NT cur_dist = get(face_L2_map, fd);

		if (max_dist < cur_dist)
		{
			if (isInsideCutpaths(sm, cutpaths, fd))
			{
				continue;
			}
			max_dist = cur_dist;
			max_distFd = fd;
		}
	}

	SM_face_descriptor sm_fd = SM_face_descriptor(max_distFd.idx());
	//randomly select a vertex in face_maxdist
	SM_halfedge_descriptor sm_hd = halfedge(sm_fd,sm);
	// the source point can not be inside the cutpaths
	// means all 
	SM_vertex_descriptor source = *CGAL::vertices_around_face(sm_hd, sm).first;
	//for (auto temp : CGAL::vertices_around_face(sm_hd, sm))
	//{
	//	bool is_inside = false;
	//	for (auto line : cutpaths)
	//	{
	//		for (auto vd : line)
	//		{
	//			if (temp == vd)
	//			{
	//				is_inside = true;
	//				break;
	//			}
	//		}
	//		if (is_inside) break;
	//	}
	//	if (!is_inside)
	//	{
	//		source = temp;
	//		break;
	//	}
	//}
	//get dijkstra map:
	SM_vertex_iterator vit, ve;
	// Associate indices to the vertices
	VertexIndexMap vertex_id_map;
	VertexIdPropertyMap vertex_index_pmap(vertex_id_map);
	std::vector<SM_vertex_descriptor> predecessor(sm.number_of_vertices());
	// and then turn it into a property map
	boost::iterator_property_map<std::vector<SM_vertex_descriptor>::iterator, VertexIdPropertyMap>
		predecessor_pmap(predecessor.begin(), vertex_index_pmap);
	std::vector<double> distance(sm.number_of_vertices());
	boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap>
		distance_pmap(distance.begin(), vertex_index_pmap);
	int index = 0;
	for (SM_vertex_descriptor vd : sm.vertices())
		vertex_id_map[vd] = index++;
	
	std::cout << "\nStart dijkstra_shortest_paths at " << source.idx() << "\n";
	boost::dijkstra_shortest_paths(
		sm, source,
		distance_map(distance_pmap)
		.predecessor_map(predecessor_pmap)
		.vertex_index_map(vertex_index_pmap));


	std::vector<SM_vertex_descriptor> path;
	SM_vertex_descriptor des_v;
	double min_distance = std::numeric_limits<double>::max();
	for (auto row : cutpaths)
	{
		for (SM_vertex_descriptor vd : row)
		{
			if (get(distance_pmap, vd) < min_distance)
			{
				min_distance = get(distance_pmap, vd);
				des_v = vd;
			}
		}
	}
	path.push_back(des_v);
	while (get(predecessor_pmap, des_v)!=source)
	{
		path.push_back(get(predecessor_pmap, des_v));
		des_v = get(predecessor_pmap, des_v);
	}
	path.push_back(source);
	cutpaths.emplace_back(path);
}



Face_NT_map get_dist_map(Seam_mesh& seam_mesh, UV_pmap& uv_map)
{
	Face_NT_map face_dist_map = get(Face_NT_tag(), seam_mesh);
	auto ppmap = get(CGAL::vertex_point, seam_mesh);
	for (face_descriptor f : seam_mesh.faces())
	{
		NT face_l2;
		{
			std::array<const Point_2*, 3> uv_points;
			std::array<const Point_3*, 3> mesh_points;

			halfedge_descriptor h = halfedge(f, seam_mesh);
			for (std::size_t i = 0; i < 3; ++i)
			{
				vertex_descriptor v = target(h, seam_mesh);

				// just to be safe in case of weird VPM returning temporaries
				boost::property_traits<UV_pmap>::reference uvp = get(uv_map, v);
				PPM_ref p = get(ppmap, v);

				uv_points[i] = &uvp;
				mesh_points[i] = &p;

				h = next(h, seam_mesh);
			}
			// Formula from Sander et al. 'Texture Mapping Progressive Meshes'
			const NT A = get_a(uv_points);
			const NT den = NT(1) / (NT(2) * A);
			const Point_3 Ss = get_ss(mesh_points, uv_points, den);
			const Point_3 St = get_st(mesh_points, uv_points, den);

			const NT a = inner_product(Ss, Ss);
			const NT c = inner_product(St, St);

			face_l2 = sqrt((a + c) / NT(2));
		}
		//      std::cout << "Face L2: " << f << " = " << compute_face_L2(f, tmesh, uvmap, ppmap) << std::endl;
		put(face_dist_map, f, face_l2);
	}
	return face_dist_map;
}


NT get_a(const std::array<const Point_2*, 3>& uv_points)
{
	NT A = (((uv_points[1]->x() - uv_points[0]->x()) * (uv_points[2]->y() - uv_points[0]->y()))
		- ((uv_points[2]->x() - uv_points[0]->x()) * (uv_points[1]->y() - uv_points[0]->y()))) / NT(2);

	CGAL_warning(A != NT(0)); // means degenerate face in the param space
	if (A == NT(0))
		return NT(1);

	return A;
}


Point_3 get_ss(const std::array<const Point_3*, 3>& mesh_points,
	const std::array<const Point_2*, 3>& uv_points,
	const NT den)
{
	const NT dt0 = uv_points[1]->y() - uv_points[2]->y();
	const NT dt1 = uv_points[2]->y() - uv_points[0]->y();
	const NT dt2 = uv_points[0]->y() - uv_points[1]->y();
	Point_3 Ss(den * (mesh_points[0]->x() * dt0 + mesh_points[1]->x() * dt1 + mesh_points[2]->x() * dt2),
		den * (mesh_points[0]->y() * dt0 + mesh_points[1]->y() * dt1 + mesh_points[2]->y() * dt2),
		den * (mesh_points[0]->z() * dt0 + mesh_points[1]->z() * dt1 + mesh_points[2]->z() * dt2));
	return Ss;
}

Point_3 get_st(const std::array<const Point_3*, 3>& mesh_points,
	const std::array<const Point_2*, 3>& uv_points,
	const NT den)
{
	const NT ds0 = uv_points[2]->x() - uv_points[1]->x();
	const NT ds1 = uv_points[0]->x() - uv_points[2]->x();
	const NT ds2 = uv_points[1]->x() - uv_points[0]->x();
	Point_3 St(den * (mesh_points[0]->x() * ds0 + mesh_points[1]->x() * ds1 + mesh_points[2]->x() * ds2),
		den * (mesh_points[0]->y() * ds0 + mesh_points[1]->y() * ds1 + mesh_points[2]->y() * ds2),
		den * (mesh_points[0]->z() * ds0 + mesh_points[1]->z() * ds1 + mesh_points[2]->z() * ds2));
	return St;
}

NT inner_product(const Point_3& pointA, const Point_3& pointB)
{
	return ((pointA.x()) * (pointB.x()) + (pointA.y()) * (pointB.y()) + (pointA.z()) * (pointB.z()));
}


NT compute_area_distortion(
	const std::vector<face_descriptor>& face_range,
	const NT A_3D,
	Seam_mesh& tmesh,
	UV_pmap&uvmap)
{


	Face_NT_map area_2DMap = get(Face_NT_tag(), tmesh);

	std::vector<NT> area_dist;
	NT A_2D = 0;

	for (face_descriptor f : face_range)
	{
		// get area in parameterised mesh
		const halfedge_descriptor h = halfedge(f, tmesh);
		const NT a_2D = abs(CGAL::area(get(uvmap, source(h, tmesh)),
			get(uvmap, target(h, tmesh)),
			get(uvmap, target(next(h, tmesh), tmesh))));
		put(area_2DMap, f, a_2D);
		A_2D += a_2D;
	}

	for (face_descriptor f : face_range)
	{
		const NT a_3D = get(m_face_areas, f);
		const NT a_2D = get(area_2DMap, f);

		area_dist.push_back(abs(a_3D / A_3D - a_2D / A_2D));
	}

	return std::accumulate(area_dist.begin(), area_dist.end(), NT(0));
}

NT initialize_faces_areas(
	const std::vector<face_descriptor>& face_range,
	Seam_mesh& tmesh)
{
	m_face_areas = get(Face_NT_tag(), tmesh);
	NT total_area = 0;

	for (face_descriptor f : face_range)
	{
		const NT f_area = PMP::face_area(f, tmesh);
		put(m_face_areas, f, f_area);
		total_area += f_area;
	}

	return total_area;
}


double get_error(Seam_mesh&seam_mesh,UV_pmap&uv_map,halfedge_descriptor&bhd)
{
	Vertex_set cc_vertices;
	std::vector<face_descriptor> cc_faces;
	cc_faces.reserve(num_faces(seam_mesh));
	SMP::internal::Containers_filler<Seam_mesh, Vertex_set> fc(seam_mesh, cc_vertices, &cc_faces);
	PMP::connected_component(face(opposite(bhd, seam_mesh), seam_mesh), seam_mesh,
		boost::make_function_output_iterator(fc));

	NT area_3D = initialize_faces_areas(cc_faces, seam_mesh);

	double error = compute_area_distortion(cc_faces, area_3D, seam_mesh, uv_map);
	return error;
}

NT compute_stretch(Seam_mesh& tmesh, Face_NT_map& face_L2_map) {
	NT stretch = 0.0;
	for (face_descriptor fd : tmesh.faces())
	{
		NT cur_dist = get(face_L2_map, fd);
		stretch += cur_dist;
	}
	return stretch;
}