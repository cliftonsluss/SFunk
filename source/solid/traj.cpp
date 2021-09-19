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
  double skin = config.j["skin"];
  std::cout << "Using skin of " << skin << std::endl;
  std::cout << "Skipping " << num_skipframes << " frames" << std::endl;
  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;
  std::string outfile = config.j["outfile"];
  int dump = config.j["dump"];


//lines below should all remove // to uncomment
// create a resultSet obect to house our results
resultSet<double> results;
//NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
// only results.avg is needed from variance00WK at thi point, that functionality
// will be broken out into a more compact method at a later date
variance00WK<double>(datafile, num_atoms, num_frames, num_skipframes, results);
// values below not needed anymore but still being reported for reference
outputfile << std::fixed << std::setprecision(7);
std::cout << "variance00= " << results.variance << std::endl;
std::cout << "std00= " << pow(results.variance,0.5) << std::endl;

size_t nbs_found;
double variance01;
std::vector<double> var_vec;

variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
  num_skipframes, num_nbs, nbs_found, variance01, outfile, skin, dump);
std::cout << "neighbor count= " << nbs_found << std::endl;
std::cout << "variance01= " << variance01 << std::endl;
std::cout << "std01= " << pow(variance01,0.5) << std::endl;
for (size_t i = 0; i < var_vec.size(); i++){
  std::cout << i+1.0 << ", " << var_vec[i] << std::endl;
}
return 0;
}
