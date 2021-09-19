#ifndef PBCPCG_H
#define PBCPCG_H
#include "structures.h"

class PBCPointCloudGenerator {
  public:
    PBCPointCloudGenerator(){};

    PBCPointCloudGenerator(simFrame<double> &frame,
        double skin, bool scaled = true, const double max_range = 10);

    PointCloud<double> GetCloud();

  private:
    PointCloud<double> cloud;
};

#endif
