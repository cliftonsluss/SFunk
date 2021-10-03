#ifndef NEIGH_LIST_GEN
#define NEIGH_LIST_GEN
#include "structures.h"
#include "PBCPointCloud_generator.h"

class NeighborListGenerator {
  public:
    NeighborListGenerator(){};

    NeighborListGenerator(simFrame<double> &avg_frame,
                          const size_t N, const int num_nbs, double skin);

    NeighborListGenerator(simFrame<double> &avg_frame,
                          const size_t N, const int num_nbs, double skin,
                          double thing);

    std::vector<std::vector<size_t>> GetListOG();
    std::vector<std::vector<size_t>> GetListSkin();
    std::vector<std::vector<size_t>> GetErrors();
    // void WriteErrors(std::string &errorfile);

  private:
    PBCPointCloudGenerator PBCPCG;
    PointCloud<double> cloud;
    std::vector<std::vector<size_t>> nbs_og;
    std::vector<std::vector<size_t>> nbs_skin;
    std::vector<std::vector<size_t>> errors;

};


#endif
