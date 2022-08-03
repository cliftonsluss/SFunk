#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <stdexcept>

#ifndef STRUCTURES_H
#define STRUCTURES_H

// common structures used throughout project

template <typename T>
struct PointCloud
{
	struct Point
	{
	  T   x,y,z;
	};

	std::vector<Point>  pts;
  std::vector<T> list;

	struct Atom {
	  int atom_num, atom_type;
	};

	std::vector<Atom>  atms;
  std::vector<int> sphere_list;
  size_t count;

	struct Bounds {
	  double min,max;
	};

	Bounds Xbounds, Ybounds, Zbounds;

	std::map<size_t,size_t> pbc_idx_map;

  // everything below is the stock point cloud from nanoflann
	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}
	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};

template<typename T>
struct Box
{
  T xmin{0}, xmax{0}, xlen{0},
    ymin{0}, ymax{0}, ylen{0},
    zmin{0}, zmax{0}, zlen{0};
};

template <typename T>
struct simFrame
{
	struct Point
	{
	  T  x,y,z;
	};

  std::vector<std::vector <double>> points;

  std::vector<size_t> idx;

	std::vector<Point>  pts;

	Box<double> box;
	struct Atom {
	  int atom_num, atom_type;
	};
	std::vector<Atom> atms;
	std::vector<double> flat;
	int num_atoms;
  int index;
};

template <typename T>
struct resultSet
{
  simFrame<T> avg;
  T variance;
  std::vector<T> var_xyz;
};

template <typename T>
struct RDF
{
  std::vector<T> r;
  std::vector<T> g_of_r;
};

// template <typename T>
// struct NeighborList
// {
//   std::vector<std::vector<T>> nbs_og;
//   std::vector<std::vector<T>> nbs_skin;
//   std::vector<std::vector<T>> errors;
// };

#endif
