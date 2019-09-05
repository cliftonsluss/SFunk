#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "test.h"

#ifndef TRAJ_READER_H
#define TRAJ_READER_H



class Trajectory {
  public:
    Trajectory(string &filename, const size_t num_atoms, const size_t header);
    void getNextFrame(simFrame<double> &frame);

  private:
    string filename;
    size_t num_atoms;
    size_t header;
    ifstream inputfile;
};

#endif
