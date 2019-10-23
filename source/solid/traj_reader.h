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
    Trajectory(std::string &filename, const size_t num_atoms, const size_t header);
    void getNextFrame(simFrame<double> &frame);
    void skipFrames(const size_t sframes);

  private:
    // string filename;
    size_t num_atoms;
    size_t header;
    std::ifstream inputfile;
};

#endif
