#include "test.h"
#include "traj_reader.h"

Trajectory::Trajectory(std::string &filename, const size_t num_atoms,
    const size_t header) {
  Trajectory::header = header;
  Trajectory::num_atoms = num_atoms;
  Trajectory::inputfile.open(filename, std::ifstream::in);
}

void Trajectory::getNextFrame(simFrame<double> &frame) {
  std::string temp;
  double x, y, z, xc, yc, zc;
  size_t n, nt;
  for (size_t i = 0; i < Trajectory::header; i++) {
    getline(inputfile, temp);
  }

  inputfile >> frame.xbox.min >> frame.xbox.max;
  double xlen = frame.xbox.max - frame.xbox.min;
  inputfile >> frame.ybox.min >> frame.ybox.max;
  double ylen = frame.ybox.max - frame.ybox.min;
  inputfile >> frame.zbox.min >> frame.zbox.max;
  double zlen = frame.zbox.max - frame.zbox.min;

  inputfile >> temp >> temp >> temp >> temp >> temp >> temp >> temp;
  frame.pts.resize(num_atoms);
  frame.atms.resize(num_atoms);
  frame.num_atoms = num_atoms;

  for (size_t i = 0; i < num_atoms; i++) {
    inputfile >> n >> nt >> xc >> yc >> zc;
    std::cout << n << std::endl;
    frame.atms[n-1].atom_num = n;
    frame.atms[n-1].atom_type = nt;
    x = xc*xlen + frame.xbox.min;
    y = yc*ylen + frame.ybox.min;
    z = zc*zlen + frame.zbox.min;
    frame.pts[n-1].x = x;
    frame.pts[n-1].y = y;
    frame.pts[n-1].z = z;
  }
  getline(inputfile, temp);
}

void Trajectory::skipFrames(const size_t sframes) {
  std::string temp;
  size_t slines = sframes*(Trajectory::header + 4 + num_atoms);
  for (size_t i = 0; i < slines; i++) {
    getline(inputfile, temp);
  }
  // statement below lets us check skipFrames functionality
  std::cout << temp << '\n';
}
