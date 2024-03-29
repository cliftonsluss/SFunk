#include "nanoflann.h"
#include "read_config.h"
#include "traj_reader.h"
#include "structures.h"
#include "PBC.h"
#include "cell2frame.h"
#include "stats.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>
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
  std::cout << "skipping " << num_skipframes << " frames" << std::endl;

  double skin = config.j["skin"];
  std::cout << "Using skin of " << skin << std::endl;

  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;

  std::string outfile = config.j["outfile"];
  std::cout << "writing running stats to " << outfile << std::endl;

  int dump = config.j["dump"];
  bool avgdump = config.j["avgdump"];
  std::string avgoutfile = config.j["avgoutfile"];

  bool nbList = config.j["nbList"];
  std::string nblfile = config.j["neighbor list file"];
  std::string avgtraj = config.j["avgtraj"];

  bool var01_only = config.j["var01_only"];
  std::string avgtrajin = config.j["avgtrajin"];





//lines below should all remove // to uncomment
// create a resultSet obect to house our results
resultSet<double> results;
simFrame<double> frame1;
std::vector<std::vector<size_t>> nbl;
//NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
// only results.avg is needed from variance00WK at thi point, that functionality
// will be broken out into a more compact method at a later date

size_t nbs_found;
double variance01;
std::vector<double> variance01_xyz;
std::vector<double> var_vec;

if (!var01_only) {
  variance00WK<double>(datafile, num_atoms,
    num_frames, num_skipframes, results);
  // values below not needed anymore but still being reported for reference
  std::cout << std::fixed << std::setprecision(16);
  // std::cout << "variance00= " << results.variance << std::endl;
  // std::cout << "std00= " << pow(results.variance,0.5) << std::endl;

  std::cout << "var00= " << results.variance << std::endl;
  std::cout << "varx00= " << results.var_xyz[0] << std::endl;
  std::cout << "vary00= " << results.var_xyz[1] << std::endl;
  std::cout << "varz00= " << results.var_xyz[2] << std::endl;
  std::cout << "std00= " << pow(results.variance,0.5) << std::endl;


  frame1 = results.avg;
  std::cout << "completed var00 calculation" << std::endl;

  if (avgdump) {
    std::cout << "writing average frame to " << avgoutfile << std::endl;
    Trajectory dump(avgoutfile);
    dump.writeFrame(frame1);
  }
} else {
  std::cout << "skipping var00 calculation" << std::endl;
  std::cout << num_atoms << std::endl;
  Trajectory traj(avgtrajin, num_atoms, 5);
  traj.getNextFrame(frame1);
  NeighborListGenerator NLG = NeighborListGenerator(
                              frame1, num_atoms, num_nbs,
                              skin);
  nbl = NLG.GetListOG();
  std::cout << "calculating var01 from nbs found in file"
  << std::endl;
  variance01kd_r_nbl<double>(datafile,
    num_atoms, num_frames,
    num_skipframes, num_nbs,
    nbs_found, variance01, variance01_xyz, outfile, skin, nbl, dump);
  std::cout << std::fixed << std::setprecision(16);
  std::cout << "variance01= " << variance01 << std::endl;
  std::cout << "var01_x= " << variance01_xyz[0] << std::endl;
  std::cout << "var01_y= " << variance01_xyz[1] << std::endl;
  std::cout << "var01_z= " << variance01_xyz[2] << std::endl;
  std::cout << "std01= " << pow(variance01,0.5) << std::endl;
  return 0;
}





// if (avgdump) {
//   Trajectory dump(avgoutfile);
//   dump.writeFrame(frame);
// }

if (nbList) {
  std::cout << "generating neighbor list for error report" <<std::endl;
  std::vector<double> pt(4);
  std::vector<double> nb(4);
  double x,y,z;
  Trajectory CellTraj(avgtraj);
  simFrame<double> frame;
  std::vector<double> len{frame1.box.xlen,
                          frame1.box.ylen,
                          frame1.box.zlen};


  std::ofstream outputfile;
  outputfile.open(nblfile, std::ofstream::out);
  NeighborListGenerator NLG = NeighborListGenerator(
                              frame1, num_atoms, num_nbs,
                              skin);
  nbl = NLG.GetListOG();
  std::vector<std::vector<size_t>> errors = NLG.GetErrors();
  std::cout << errors.size() << std::endl;
  if (errors.size() > 0){
    outputfile << errors.size() <<
      " pairs with mismatched neighbors" << std::endl;
    for (int i = 0; i < errors.size(); i++){
      outputfile << errors[i][0]+1 << ", "
        << errors[i][1]+1 << "\n" << std::endl;
    }
    outputfile << "\n" << std::endl;
    for (int i = 0; i < errors.size(); i++){

      x = frame1.pts[errors[i][0]].x;
      y = frame1.pts[errors[i][0]].y;
      z = frame1.pts[errors[i][0]].z;
      pt = {errors[i][0]+1.0,x,y,z};
      PBC pbc1(&pt[1],len);
      std::vector<std::vector<double>> nbs1;

      outputfile << "neighbors of " << errors[i][0]+1 << std::endl;

      for (int j = 0; j < num_nbs; j++){
        x = frame1.pts[nbl[errors[i][0]][j]].x;
        y = frame1.pts[nbl[errors[i][0]][j]].y;
        z = frame1.pts[nbl[errors[i][0]][j]].z;
        nb = {nbl[errors[i][0]][j]+1.0,x,y,z};
        pbc1.minimum_image(&nb[1]);
        nbs1.push_back(nb);

        outputfile << "ParticleIdentifier==" << nbl[errors[i][0]][j]+1
        << " ||" << std::endl;
      }
      outputfile << "\n";
      frame = cell2frame(pt, nbs1);
      CellTraj.writeFrame(frame);

      x = frame1.pts[errors[i][1]].x;
      y = frame1.pts[errors[i][1]].y;
      z = frame1.pts[errors[i][1]].z;

      pt = {errors[i][1]+1.0,x,y,z};
      PBC pbc2(&pt[1],len);
      std::vector<std::vector<double>> nbs2;

      outputfile << "ParticleIdentifier==" << errors[i][0]+1 << "\n" << std::endl;
      outputfile << "******************************************" << std::endl;

      outputfile << "neighbors of " << errors[i][1]+1 << std::endl;

      for (int j = 0; j < num_nbs; j++){
        x = frame1.pts[nbl[errors[i][1]][j]].x;
        y = frame1.pts[nbl[errors[i][1]][j]].y;
        z = frame1.pts[nbl[errors[i][1]][j]].z;
        nb = {nbl[errors[i][1]][j]+1.0,x,y,z};
        pbc2.minimum_image(&nb[1]);
        nbs2.push_back(nb);

        outputfile << "ParticleIdentifier==" << nbl[errors[i][1]][j]+1
        << " ||" << std::endl;
      }
      outputfile << "\n";
      frame = cell2frame(pt, nbs2);
      CellTraj.writeFrame(frame);


      outputfile << "ParticleIdentifier==" << errors[i][1]+1 << "\n" << std::endl;
      outputfile << "******************************************" << std::endl;
    }
  }


}

// if (var01_only) {
//   std::cout << "calculating var01 from nbs found in file" << std::endl;
//   variance01kd_r_nbl<double>(datafile, num_atoms, num_frames,
//     num_skipframes, num_nbs, nbs_found, variance01, outfile, skin, nbl, dump);
//     std::cout << "variance01= " << variance01 << std::endl;
//     std::cout << "std01= " << pow(variance01,0.5) << std::endl;
// } else {
std::cout << "calculating var01 from var00 output" << std::endl;
variance01kd_r<double>(datafile, frame1, num_atoms, num_frames,
  num_skipframes, num_nbs, nbs_found, variance01, variance01_xyz,
  outfile, skin, dump);
// std::cout << "neighbor count= " << nbs_found << std::endl;
std::cout << std::fixed << std::setprecision(16);
std::cout << "variance01= " << variance01 << std::endl;
std::cout << "var01_x= " << variance01_xyz[0] << std::endl;
std::cout << "var01_y= " << variance01_xyz[1] << std::endl;
std::cout << "var01_z= " << variance01_xyz[2] << std::endl;
std::cout << "std01= " << pow(variance01,0.5) << std::endl;
// }

// for (size_t i = 0; i < var_vec.size(); i++){
//   std::cout << i+1.0 << ", " << var_vec[i] << std::endl;
// }
return 0;
}
