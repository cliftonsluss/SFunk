#include "test.h"


Trajectory::Trajectory(string &filename, const size_t num_atoms, const size_t header) {
  Trajectory::header = header;
  Trajectory::num_atoms = num_atoms;
  Trajectory::inputfile.open(filename, std::ifstream::in);
}

void Trajectory::getNextFrame(simFrame<double> &frame) {
  string temp;
  double x, y, z, xc, yc, zc;
  int n, nt;
  for (int i = 0; i < Trajectory::header; i++) {
    getline(inputfile, temp);
  }

  inputfile >> frame.xbox.min >> frame.xbox.max;
  double xlen = frame.xbox.max - frame.xbox.min;
  //cout << xlen << endl;
  inputfile >> frame.ybox.min >> frame.ybox.max;
  double ylen = frame.ybox.max - frame.ybox.min;
  inputfile >> frame.zbox.min >> frame.zbox.max;
  double zlen = frame.zbox.max - frame.zbox.min;

  inputfile >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
  frame.pts.resize(num_atoms);
  frame.atms.resize(num_atoms);
  frame.num_atoms = num_atoms;
  
  for (int i = 0; i < num_atoms; i++) {
    inputfile >> n >> nt >> xc >> yc >> zc;
    //cout << n << endl;
    frame.atms[n-1].atom_num = n;
    frame.atms[n-1].atom_type = nt;
    x = xc*xlen + frame.xbox.min;
    //cout << x << endl;
    y = yc*ylen + frame.ybox.min;
    z = zc*zlen + frame.zbox.min;
    frame.pts[n-1].x = x;
    frame.pts[n-1].y = y;
    frame.pts[n-1].z = z;
  }
  getline(inputfile, temp);
}

