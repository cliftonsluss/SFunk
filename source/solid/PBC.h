#include <vector>
#include "structures.h"

#ifndef PBC_H
#define PBC_H

class PBC {
public:
  PBC(double *pt, std::vector<double> &len);
  void minimum_image(double *ipt);
private:
  double *pt;
  std::vector<double> len; //{box.xlen, box.ylen, box.zlen};
};

#endif
