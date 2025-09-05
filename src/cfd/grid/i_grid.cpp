#include "i_grid.hpp"

using namespace cfd;

///////////////////////////////////////////////////////////////////////////////
// Cache
///////////////////////////////////////////////////////////////////////////////
void IGrid::Cache::clear(){
	boundary_faces.clear();
	boundary_points.clear();
}

void IGrid::Cache::need_boundary_faces(const IGrid& grid){
	if (boundary_faces.size() > 0){
		return;
	}
	for (size_t iface=0; iface<grid.n_faces(); ++iface){
		std::array<size_t, 2> cc = grid.tab_face_cell(iface);
		if (cc[0] == INVALID_INDEX || cc[1] == INVALID_INDEX){
			boundary_faces.push_back(iface);
		}
	}
}

void IGrid::Cache::need_boundary_points(const IGrid& grid){
	if (boundary_points.size() > 0){
		return;
	}
	need_boundary_faces(grid);
	std::set<size_t> points;
	for (size_t iface: boundary_faces){
		for (size_t ipoint: grid.tab_face_point(iface)){
			points.insert(ipoint);
		}
	}
	boundary_points = std::vector<size_t>(points.begin(), points.end());
}


///////////////////////////////////////////////////////////////////////////////
// IGrid
///////////////////////////////////////////////////////////////////////////////
std::vector<size_t> IGrid::boundary_faces() const{
	_cache.need_boundary_faces(*this);
	return _cache.boundary_faces;
}

std::vector<size_t> IGrid::boundary_points() const{
	_cache.need_boundary_points(*this);
	return _cache.boundary_points;
}
