#include "structures.h"

simFrame<double> cell2frame(std::vector<double> pt,
                            std::vector<std::vector<double>> &nbs){
  simFrame<double> frame;
  frame.pts.resize(nbs.size()+1);
  frame.idx.resize(nbs.size()+1);
  // frame.pts[0].x = pt[1];
  // frame.pts[0].y = pt[2];
  // frame.pts[0].z = pt[3];
  frame.pts[0].x = 0.0;
  frame.pts[0].y = 0.0;
  frame.pts[0].z = 0.0;
  frame.idx[0] = pt[0];

  for (int i = 1; i < nbs.size()+1; i++){
    frame.idx[i] = nbs[i-1][0];
    // frame.pts[i].x = nbs[i-1][1];
    // frame.pts[i].y = nbs[i-1][2];
    // frame.pts[i].z = nbs[i-1][3];
    frame.pts[i].x = nbs[i-1][1] - pt[1];
    frame.pts[i].y = nbs[i-1][2] - pt[2];
    frame.pts[i].z = nbs[i-1][3] - pt[3];
    frame.box.xmax = (frame.pts[i].x > frame.box.xmax)
          ? frame.pts[i].x:frame.box.xmax;
    frame.box.xmin = (frame.pts[i].x < frame.box.xmin)
          ? frame.pts[i].x:frame.box.xmin;
    frame.box.ymax = (frame.pts[i].y > frame.box.ymax)
          ? frame.pts[i].y:frame.box.ymax;
    frame.box.ymin = (frame.pts[i].y < frame.box.ymin)
          ? frame.pts[i].y:frame.box.ymin;
    frame.box.zmax = (frame.pts[i].z > frame.box.zmax)
          ? frame.pts[i].z:frame.box.zmax;
    frame.box.zmin = (frame.pts[i].z < frame.box.zmin)
          ? frame.pts[i].z:frame.box.zmin;
  }
  frame.box.xlen = frame.box.xmax - frame.box.xmin;
  frame.box.ylen = frame.box.ymax - frame.box.ymin;
  frame.box.zlen = frame.box.zmax - frame.box.zmin;
  frame.num_atoms = nbs.size() +1;

  return frame;

};
