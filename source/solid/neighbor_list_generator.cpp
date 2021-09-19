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
      // for any given neighbor index (neigh_idx) that is greater than the
      // current atom's index (idx) we want to add it to the list of neighbor
      // indexs (neigh_idxs) the greater than condition ensures we don't include
      // duplicate pairs
      // if (neigh_idx > idx){
      //   if (neigh_idx < N){
      //     NList.idx.push_back(idx);
      //     NList.nbs.push_back(neigh_idx);
      //     // std::cout << idx << " has neighbor " << neigh_idx <<std::endl;
      //     NList.nnbs++;
      //   }
      //   else {
      //     if (cloud.pbc_idx_map[neigh_idx] > idx){
      //       NList.idx.push_back(idx);
      //       NList.nbs.push_back(cloud.pbc_idx_map[neigh_idx]);
      //       // std::cout << idx << " has skin neighbor " << neigh_idx <<std::endl;
      //       NList.nnbs++;
      //     }
      //   }
      // }
      // if (idx==4){
      //   std::cout << idx << " has neighbor " << neigh_idx <<std::endl;
      // }

      nbs_skin[idx].push_back(neigh_idx);

      if (neigh_idx < N){
        // NList.idx.push_back(idx);
        // NList.nbs.push_back(neigh_idx);
        //
        // NList.nnbs++;
        nbs_og[idx][k] = neigh_idx;
      }
      else {
        // if (cloud.pbc_idx_map[neigh_idx] > idx){
          // NList.idx.push_back(idx);
          // NList.nbs.push_back(cloud.pbc_idx_map[neigh_idx]);
          // // std::cout << idx << " has skin neighbor " << neigh_idx <<
          // // " which has original index " << cloud.pbc_idx_map[neigh_idx] <<std::endl;
          // NList.nnbs++;
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





      // if (neigh_idx > idx) {
      //   // is the neighbor index found in the skin map
      //   if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
      //     if (cloud.pbc_idx_map[neigh_idx] > idx) {
      //       NList.idx.push_back(idx);
      //       NList.nbs.push_back(cloud.pbc_idx_map[neigh_idx]);
      //       //std::cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << std::endl;
      //       //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
      //       NList.nnbs++;
      //     }
      //   }
      //   // if the index is not found in the skin map the pair
      //   // can be stored as is
      //   else {
      //     NList.idx.push_back(idx);
      //     NList.nbs.push_back(neigh_idx);
      //     //std::cout << idx << ", " << neigh_idx << std::endl;
      //     //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
      //     NList.nnbs++;
      //   }
      // }




    }
  }
}

// std::cout << "atom 4 at x = " << cloud.pts[4].x << " y = " << cloud.pts[4].y <<
// " = z " << cloud.pts[4].z << std::endl;
//
// std::cout << "atom 2596 at x = " << cloud.pts[2596].x << " y = " << cloud.pts[2596].y <<
// " = z " << cloud.pts[2596].z << std::endl;
// int a,b,a1,b1,k,m;
//
//
// for (int i =0; i < N*num_nbs; i++){
//   a = NList.idx[i];
//   b = NList.nbs[i];
//   k = 0;
//   for (int j =0; j < N*num_nbs; j++){
//     a1 = NList.idx[j];
//     b1 = NList.nbs[j];
//     if (a1 == b){
//       if (b1 == a){
//       k = 1;
//       }
//     }
//   }
//
//   if (k == 0){
//     std::cout << "a = " << a << ", b = " << b << std::endl;
//   }
//
// }

};

// NeighborListGenerator::NeighborListGenerator(simFrame<double> &avg_frame,
//                 const size_t N,
//                 const int num_nbs,
//                 double skin,
//                 double thing){
//   PointCloud<double> cloud;
//   populatePointCloudSkin(cloud, avg_frame, skin);
//
//   // set up KDtree adaptor
//   typedef KDTreeSingleIndexAdaptor<
//    NN_Adaptor<double, PointCloud<double> >,
//    PointCloud<double>,
//    3
//    > my_kd_tree_t;
//
//   // create an index for our adaptor and build index
//   my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(12) );
//   index.buildIndex();
//   size_t num_results = num_nbs;
//   int idx;
//   int neigh_idx;
//   // NeighborList<size_t> NList {std::vector<int>{},std::vector<int>{},0};
//
//   // nnbs = 0;
//
//   // we will use each atom inside the threshold as a query point to find
//   // the nearest neighbors
//   for (int j = 0; j < cloud.list.size(); j++)
//   {
//     idx = cloud.list[j];
//     // std::cout << idx << "\n";
//     const double query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
//         cloud.pts[idx].z };
//     {
//     std::vector<size_t> ret_index(num_results);
//     std::vector<double> out_dist(num_results);
//
//     // pass reference to location of query_pt
//     num_results = index.knnSearch(&query_pt[0], num_nbs, &ret_index[0],
//         &out_dist[0]);
//     // std::cout << "num_results= " << num_results << std::endl;
//     // std::cout << "num_nbs= " << num_nbs << std::endl;
//
//     // catch any wrong number of neighbors produced by nanoflann
//     if (num_results != num_nbs) {
//       try {
//         throw num_results;
//       }
//       catch (size_t i ) {
//         std::cout << "Found " << i << " nearest neighbors, expected "
//         << num_nbs << std::endl;
//       }
//     }
//     ret_index.resize(num_results);
//     out_dist.resize(num_results);
//     int temp;
//     for (size_t k = 0; k < num_results; k++) {
//       neigh_idx = ret_index[k];
//       // std::cout << neigh_idx << "\n";
//       if (neigh_idx > idx){
//         NList.idx.push_back(idx);
//         NList.nbs.push_back(neigh_idx);
//         // NList.nnbs++;
//       }
//     }
//     }
//   }
// };

std::vector<std::vector<size_t>> NeighborListGenerator::GetListOG(){
  return nbs_og;
};

std::vector<std::vector<size_t>> NeighborListGenerator::GetListSkin(){
  return nbs_skin;
};

std::vector<std::vector<size_t>> NeighborListGenerator::GetErrors(){
  return errors;
};
