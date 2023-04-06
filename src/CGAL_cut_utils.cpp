#include "CGAL_cut_utils.h"


bool isInsideCutpaths(
	const Surface_mesh&sm,
	std::vector<std::vector<SM_vertex_descriptor>>& cutpaths,
	face_descriptor fd)
{
	SM_face_descriptor sm_fd = SM_face_descriptor(fd.id());
	//randomly select a vertex in face_maxdist
	SM_halfedge_descriptor sm_hd = halfedge(sm_fd, sm);
	bool is_inside = false;
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
					is_inside = true;
					break;
				}
			}
			if (is_inside) break;
		}
		if (!is_inside)
		{
			break;
		}
	}
	return is_inside;
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

	SM_face_descriptor sm_fd = SM_face_descriptor(max_distFd.id());
	//randomly select a vertex in face_maxdist
	SM_halfedge_descriptor sm_hd = halfedge(sm_fd,sm);
	// the source point can not be inside the cutpaths
	// means all 
	SM_vertex_descriptor source;
	for (auto temp : CGAL::vertices_around_face(sm_hd, sm))
	{
		bool is_inside = false;
		for (auto line : cutpaths)
		{
			for (auto vd : line)
			{
				if (temp == vd)
				{
					is_inside = true;
					break;
				}
			}
			if (is_inside) break;
		}
		if (!is_inside)
		{
			source = temp;
			break;
		}
	}
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
	int min_distance = INT_MAX;
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




