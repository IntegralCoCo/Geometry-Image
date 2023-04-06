#include "cgal_IO.h"

bool CGAL_read_OBJ(std::istream& is, Surface_mesh& sm)
{
	std::string item, line, col;
	Point_3 p;
	Vector_3 v;

	std::vector<Vertex> vertexmap;
	std::vector<Vector_3> normals;
	std::vector<Vertex> vr(3);

	VProp_geom vpm = sm.points();
	auto vnormal = sm.add_property_map<Vertex, Vector_3>("v:normal", Vector_3(0, 0, 0)).first;
	auto vcolor = sm.add_property_map<Vertex, CGAL::Color>("v:color", CGAL::Color(0, 0, 0)).first;

	//is >> CGAL::sm_skip_comments;
	size_t vi;
	while (getline(is, line)) {
		std::istringstream iss(line);
		iss >> item;
		if (item == "#")
		{
			continue;
		}
		if (item == "vn")
		{
			iss >> v;
			normals.push_back(v);
		}
		else if (item == "v")
		{
			iss >> p;
			Vertex vi = sm.add_vertex();
			vpm[vi] = p;
			vertexmap.push_back(vi);
			// For color
			std::getline(iss, col);
			std::istringstream isss(col);
			if (isss >> item)
			{
				vcolor[vi] = CGAL::File_scanner_OFF::get_color_from_line(std::istringstream(col));
			}
		}
		else if (item == "f") {
			vr.resize(3);
			for (std::size_t j = 0; j < 3; j++) {
				iss >> item;
				std::istringstream isss(item.substr(0, item.find("\\")));
				isss >> vi;
				vr[j] = vertexmap[vi - 1];
			}
			Face fi = sm.add_face(vr);
			if (fi == sm.null_face())
			{
				is.setstate(std::ios::failbit);
				sm.clear();
				return false;
			}
		}
		else
			continue;
	}
	if (!normals.empty())
	{
		if (normals.size() == vertexmap.size())
		{
			for (int i = 0; i < normals.size(); ++i)
				vnormal[vertexmap[i]] = normals[i];
		}
		else
		{
			is.setstate(std::ios::failbit);
			return false;
		}
	}
	return is.good();
}

bool outfileExists(fs::path outFilePath, const int size, std::string printDesc) {
	if (fs::exists(outFilePath)) {
		if (fs::file_size(outFilePath) > size) {
			std::cout << printDesc << std::flush;
			return true;
		}
		else
			return false;
	}
	else
		return false;
}

bool infileExists(fs::path inFilePath, const int size, std::string errDesc) {
	if (fs::exists(inFilePath)) {
		//check for size of file
		if (fs::file_size(inFilePath) <= size) {
			std::cerr << "\t Empty" << errDesc << "\n";
			return false;
		}
		else
			return true;
	}
	else {
		std::cerr << "\t No " << errDesc << "\n";
		return false;
	}
}

bool meshLoader(fs::path meshFile, Surface_mesh& loadedMesh, std::string fileDesc, bool bdebug) {
	// check if the file exists
	if (fs::exists(meshFile) && fs::is_regular_file(meshFile)) {
		fs::ifstream in_fs(meshFile);
		if (!in_fs) {
			std::cerr << "\t Unable to create fs for " << fileDesc << std::endl;
			return false;
		}
		//CGAL_read_OBJ(in_fs, loadedMesh);
		try {
			if (!(in_fs >> loadedMesh)) {
				std::cerr << std::setw(20) << "\t Unable to read " << fileDesc << std::endl;
				return false;
			}
		}
		catch (...) {
			std::cerr << std::setw(20) << "\t Unable to read " << fileDesc << std::endl;
			return false;
		}

		if (loadedMesh.number_of_vertices() == 0) {
			std::cerr << "\t" << fileDesc << "has no vertices" << std::endl;
			return false;
		}
		else {
			if (bdebug)
				std::cout << "Loaded Mesh " << meshFile << " has " << loadedMesh.number_of_vertices() << " Vertices "
				<< std::endl;
		}
		in_fs.close();
		return true;
	}
	else {
		std::cerr << "\t" << fileDesc << "doesn't exists\n";
		return false;
	}
}




