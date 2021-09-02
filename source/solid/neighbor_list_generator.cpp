#include "test.h"
#include "neighbor_list_generator.h"
#include <vector>
#include <stdexcept>
#include <sstream>


NeighborListGenerator::NeighborListGenerator(simFrame<double> &avg_frame,
                const size_t N, const int num_nbs, double skin){

  // Create point cloud and populate with average frame utilizing
  // periodic boundary conditions
  PointCloud<double> cloud;
  populatePointCloudPBC(cloud, avg_frame, skin);

  // set up KDtree adaptor
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<double, PointCloud<double> >,
   PointCloud<double>,
   3
   > my_kd_tree_t;

  // create an index for our adaptor and build index
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  index.buildIndex();
  size_t num_results = num_nbs;
  int idx;
  int neigh_idx;
  // NeighborList<size_t> NList {std::vector<int>{},std::vector<int>{},0};

  NList.nnbs = 0;

  // we will use each atom as a query point to find nearest neighbors
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const double query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<double> out_dist(num_results);

    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0],
        &out_dist[0]);
    if (num_results != num_nbs) {
      try {
        throw num_results;
      }
      catch (size_t i ) {
        std::cout << "Found " << i << "neighbors, expected "
        << num_nbs << std::endl;
      }
    }
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    // downselect neighbors and restore indexes from skin map
    // count number of actual neighbors
    for (size_t k = 0; k < num_results; k++) {
      neigh_idx = ret_index[k];
      // for any given neighbor index (neigh_idx) that is greater than the
      // current atom's index (idx) we want to add it to the list of neighbor
      // indexs (neigh_idxs) the greater than condition ensures we don't include
      // duplicate pairs
      if (neigh_idx > idx) {
        // is the neighbor index found in the skin map
        if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
          if (cloud.pbc_idx_map[neigh_idx] > idx) {
            NList.idx.push_back(idx);
            NList.nbs.push_back(cloud.pbc_idx_map[neigh_idx]);
            //std::cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << std::endl;
            //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
            NList.nnbs++;
          }
        }
        // if the index is not found in the skin map the pair
        // can be stored as is
        else {
          NList.idx.push_back(idx);
          NList.nbs.push_back(neigh_idx);
          //std::cout << idx << ", " << neigh_idx << std::endl;
          //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
          NList.nnbs++;
        }
      }
    }
  }
}
};

NeighborList<size_t> NeighborListGenerator::GetList(){
  return NList;
};
