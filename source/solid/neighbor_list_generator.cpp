#include "test.h"
#include "neighbor_list_generator.h"
#include <vector>
#include <stdexcept>
#include <sstream>


NeighborListGenerator::NeighborListGenerator(simFrame<double> &avg_frame,
                const size_t N, const int num_nbs, double skin){

  // Create point cloud and populate with average frame utilizing
  // periodic boundary conditions
  // PointCloud<double> cloud;
  // populatePointCloudPBC(cloud, avg_frame, skin);
  PBCPCG = PBCPointCloudGenerator(avg_frame, skin);
  cloud = PBCPCG.GetCloud();

  // set up KDtree adaptor
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<double, PointCloud<double> >,
   PointCloud<double>,
   3
   > my_kd_tree_t;

  // create an index for our adaptor and build index
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(6) );
  index.buildIndex();
  size_t num_results = num_nbs;
  size_t idx;
  size_t neigh_idx;

  nbs_og.resize(N, std::vector<size_t>(num_nbs));
  nbs_skin.resize(N, std::vector<size_t>(num_nbs));
  std::vector<size_t> buff(2);

  // we will use each atom as a query point to find nearest neighbors
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const double query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<double> out_dist(num_results);

    // pass reference to location of query_pt
    num_results = index.knnSearch(&query_pt[0], num_nbs, &ret_index[0],
        &out_dist[0]);
    if (num_results != num_nbs) {
      try {
        throw num_results;
      }
      catch (size_t i ) {
        std::cout << "Found " << i << " nearest neighbors, expected "
        << num_nbs << std::endl;
      }
    }
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    // downselect neighbors and restore indexes from skin map
    // count number of actual neighbors
    for (size_t k = 0; k < num_results; k++) {
      neigh_idx = ret_index[k];
      nbs_skin[idx].push_back(neigh_idx);

      if (neigh_idx < N){
        nbs_og[idx][k] = neigh_idx;
      }
      else {
        nbs_og[idx][k] = cloud.pbc_idx_map[neigh_idx];
      }
      if (neigh_idx < idx){
        if(find (nbs_og[neigh_idx].begin(),
                 nbs_og[neigh_idx].end(),
                 idx) == nbs_og[neigh_idx].end()){
          buff[0] = neigh_idx;
          buff[1] = idx;
          errors.push_back(buff);
        }
      }
    }
  }
}
};

// void NeighborListGenerator::WriteErrors(std::string &errorfile) {
//
//
// }

std::vector<std::vector<size_t>> NeighborListGenerator::GetListOG(){
  return nbs_og;
};

std::vector<std::vector<size_t>> NeighborListGenerator::GetListSkin(){
  return nbs_skin;
};

std::vector<std::vector<size_t>> NeighborListGenerator::GetErrors(){
  return errors;
};
