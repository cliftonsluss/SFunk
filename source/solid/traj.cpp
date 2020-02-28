#include "nanoflann.h"
#include "read_config.h"
#include "traj_reader.h"
#include "test.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>

int main(int argc, char *argv[]) {
  // read parameters from json file

  std::string config_file = argv[1];
  std::cout << "Reading configuration file..." << std::endl;
  Read_config config(config_file);
  config.get_config();
  int num_nbs = config.j["neighbors"];
  std::cout << num_nbs << " nearest neighbors" << std::endl;
  int num_atoms = config.j["atoms"];
  std::cout << num_atoms << " atoms" << std::endl;
  int num_frames = config.j["frames"];
  std::cout << num_frames << " frames" << std::endl;
  int num_skipframes = config.j["skipframes"];
  std::cout << "Skipping " << num_skipframes << " frames" << std::endl;
  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;


//lines below should all remove // to uncomment
// create a resultSet obect to house our results
resultSet<double> results;
//NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
// only results.avg is needed from variance00WK at thi point, that functionality
// will be broken out into a more compact method at a later date
variance00WK<double>(datafile, num_atoms, num_frames, num_skipframes, results);
// values below not needed anymore but still being reported for reference
std::cout << "variance00= " << results.variance << std::endl;
std::cout << "std00= " << pow(results.variance,0.5) << std::endl;

size_t nbs_found;
float variance01;

variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
  num_skipframes, num_nbs, nbs_found, variance01);
std::cout << "neighbor count= " << nbs_found << std::endl;
std::cout << "variance01= " << variance01 << std::endl;
std::cout << "std01= " << pow(variance01,0.5) << std::endl;
return 0;
}
