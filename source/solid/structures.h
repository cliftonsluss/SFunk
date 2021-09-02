#ifndef STRUCTURES_H
#define STRUCTURES_H

// common structures used throughout project

template <typename T>
struct simFrame
{
	struct Point
	{
	  T  x,y,z;
	};
	std::vector<Point>  pts;
	struct Box {
	  T min, max;
	};
	Box xbox, ybox, zbox;
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
};

#endif
