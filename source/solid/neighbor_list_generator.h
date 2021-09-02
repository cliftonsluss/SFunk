#ifndef NEIGH_LIST_GEN
#define NEIGH_LIST_GEN
#include "structures.h"

template <typename T>
struct NeighborList
{
  std::vector<size_t> idx;
  std::vector<size_t> nbs;
  T nnbs;
};



class NeighborListGenerator {
  public:
    NeighborListGenerator(){};

    NeighborListGenerator(simFrame<double> &avg_frame,
                          const size_t N, const int num_nbs, double skin);

    NeighborList<size_t> GetList();

  private:
    NeighborList<size_t> NList;

};


#endif
