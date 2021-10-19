#include "PBC.h"
#include "math.h"


PBC::PBC(double *pt, std::vector<double> &len){
  PBC::pt = pt;
  PBC::len = len;
};

void PBC::minimum_image(double *ipt){
  // int x = (ipt.size() > 3)?1:0;
  for (int i = 0;i < 3;i++){
    if ((pt[i]-ipt[i]) < (-len[i]*0.5)){
      ipt[i] = ipt[i] - len[i];
    }
    if ((pt[i]-ipt[i]) > (len[i]*0.5)){
      ipt[i] = ipt[i] + len[i];
    }
  }
  // std::cout << ipt[0] << "," << ipt[1] << "," << ipt[2] << std::endl;
}

double PBC::minimum_image_L2_distance(double *ipt){
  double distance = 0.0;
  double temp = 0.0;
  for (int i = 0;i < 3;i++){
    if ((pt[i]-ipt[i]) < (-len[i]*0.5)){
      temp = pt[i] - (ipt[i] - len[i]);
    }
    else if ((pt[i]-ipt[i]) > (len[i]*0.5)){
      temp = pt[i] - (ipt[i] + len[i]);
    } else {
      temp = pt[i] - ipt[i];
    }
    distance = distance + pow(temp,2);
  }
  return pow(distance, 0.5);
}

std::vector<double> PBC::minimum_image_xyz_distance(double *ipt){
  std::vector<double> distance{0,0,0};
  double temp = 0;
  for (int i = 0;i < 3;i++){
    if ((pt[i]-ipt[i]) < (-len[i]*0.5)){
      temp = abs(pt[i] - (ipt[i] - len[i]));
    }
    else if ((pt[i]-ipt[i]) > (len[i]*0.5)){
      temp = abs(pt[i] - (ipt[i] + len[i]));
    } else {
      temp = abs(pt[i] - ipt[i]);
    }
    std::cout << temp << std::endl;
    distance[i] = temp;
  }
  return distance;
}
