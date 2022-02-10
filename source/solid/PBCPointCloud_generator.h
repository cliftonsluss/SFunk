#ifndef PBCPCG_H
#define PBCPCG_H
#include "structures.h"

class PBCPointCloudGenerator {
  public:
    PBCPointCloudGenerator(){};

    PBCPointCloudGenerator(simFrame<double> &frame,
        double skin, bool scaled = true);

    PointCloud<double> GetCloud();

  private:
    PointCloud<double> cloud;
};

#endif
