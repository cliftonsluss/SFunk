#include "PBC.h"


PBC::PBC(double *pt, std::vector<double> &len){
  PBC::pt = pt;
  PBC::len = len;
};

void PBC::minimum_image(double *ipt){
  // int x = (ipt.size() > 3)?1:0;
  for (int i = 0;i < 3;i++){
    if ((ipt[i]-pt[i]) < (-len[i]*0.5)){
      ipt[i] = ipt[i] + len[i];
    }
    if ((ipt[i]-pt[i]) > (len[i]*0.5)){
      ipt[i] = ipt[i] - len[i];
    }
  }
}
