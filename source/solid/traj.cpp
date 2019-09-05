#include "nanoflann.h"
#include "test.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace nanoflann;
// function to calculate mean and standard deviation given a set of data in the
// std::vector "data" returns a std::vector {mean, standard deviation}
std::vector<double> vectorStats(std::vector<double> data) {
  double sum = accumulate(data.begin(), data.end(), 0.0);
  double mean = sum/data.size();
  std::vector<double> diff(data.size());
  transform(data.begin(), data.end(), diff.begin(),
      bind2nd(std::minus<double>(), mean));
  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = sqrt(sq_sum / data.size());
  std::vector<double> out = {mean, stdev};
  return out;
}
template <typename num_t>
void NNeighbors(std::string &filename, const size_t N, const int num_frames,
    const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 2.5;
  size_t header = 5;
  int n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist,
      zdist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame;
  for (int i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    populatePointCloudPBC(cloud, frame, skin, true);
    typedef KDTreeSingleIndexAdaptor<
      //L2_Simple_Adaptor<num_t, PointCloud<num_t> >,
      NN_Adaptor<num_t, PointCloud<num_t> >,
      PointCloud<num_t>,
      3 //dimensionality of data
      > my_kd_tree_t;
    my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
    index.buildIndex();
    size_t num_results = num_nbs;
    size_t idx;
    size_t neigh_idx;
    size_t nbs = 0;
    std::vector<size_t> idxs;
    std::vector<size_t> neigh_idxs;
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    for (int j = 0; j < N; j++)
    {
      idx = j;
      const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
          cloud.pts[idx].z };
      {
      std::vector<size_t> ret_index(num_results);
      std::vector<num_t> out_dist(num_results);
      num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0],
          &out_dist[0]);
      ret_index.resize(num_results);
      out_dist.resize(num_results);
      // downselect neighbors and restore indexes from skin map
      // count number of actual neighbors
      for (size_t k = 0; k < num_results; k++) {
        neigh_idx = ret_index[k];
        if (neigh_idx > idx) {
          // is the neighbor index found in the skin map
          if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
  	        if (cloud.pbc_idx_map[neigh_idx] > idx) {
  	          idxs.push_back(idx);
  	          neigh_idxs.push_back(cloud.pbc_idx_map[neigh_idx]);
              //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
              nbs++;
  	        }
  	      }
          // if the index is not found in the skin map the pair
          // can be stored as is
  	      else {
  	        idxs.push_back(idx);
  	        neigh_idxs.push_back(neigh_idx);
            //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
            nbs++;
  	      }
        }
      }
      }
      for (size_t k = 0; k < nbs; k++) {
        xa = frame.pts[idxs[k]].x;
        xb = frame.pts[neigh_idxs[k]].x;
        ya = frame.pts[idxs[k]].y;
        yb = frame.pts[neigh_idxs[k]].y;
        za = frame.pts[idxs[k]].z;
        zb = frame.pts[neigh_idxs[k]].z;
        if ((xa-xb) < (-xlen*0.5)) {
  	      xb = xb - xlen;
        }
        if ((xa-xb) > (xlen*0.5)){
  	      xb = xb + xlen;
        }
        if ((ya-yb) < (-ylen*0.5)){
  	      yb = yb - ylen;
        }
        if ((ya-yb) > (ylen*0.5)){
  	      yb = yb + ylen;
        }
        if ((za-zb) < (-zlen*0.5)){
  	      zb = zb - zlen;
        }
        if ((za-zb) > (zlen*0.5)){
  	      zb = zb + zlen;
        }
        xdist = std::abs(xa-xb);
        ydist = std::abs(ya-yb);
        zdist = std::abs(za-zb);
        old_avg = avg;
        avg = old_avg + (xdist - old_avg)/n;
        diff_sqrd = diff_sqrd + ((xdist - old_avg)*(xdist - avg));
        n++;
        old_avg = avg;
        avg = old_avg + (ydist - avg)/n;
        diff_sqrd = diff_sqrd + ((ydist - old_avg)*(ydist - avg));
        n++;
        old_avg = avg;
        avg = old_avg + (zdist - avg)/n;
        diff_sqrd = diff_sqrd + ((zdist - old_avg)*(zdist - avg));
        n++;
      }
      idxs.clear();
      neigh_idxs.clear();
      nbs = 0;
    }
  }
  std::cout << "average= " << avg << std::endl;
  variance = diff_sqrd/(n-1);
  std::cout << "n pairs= " << n << std::endl;
  std::cout << "variance01= " << variance << std::endl;
  std::cout << "std01= " << pow(variance, 0.5) << std::endl;
}

template <typename num_t>
void variance01kd(std::string &filename, simFrame<num_t> &avg_frame, const size_t N,
    const int num_frames, const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 2.5;
  size_t header = 5;
  int n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist, zdist,
      variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame;
  //std::cout << "populating point cloud" << std::endl;
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  //std::cout << "building point cloud" << std::endl;
  index.buildIndex();
  size_t num_results = num_nbs;
  int idx;
  int neigh_idx;
  int nbs = 0;
  std::vector<double> idxs;
  std::vector<double> neigh_idxs;
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<num_t> out_dist(num_results);
    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0],
        &out_dist[0]);
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    // downselect neighbors and restore indexes from skin map
    // count number of actual neighbors
    for (size_t k = 0; k < num_results; k++) {
      neigh_idx = ret_index[k];
      if (neigh_idx > idx) {
        // is the neighbor index found in the skin map
	      if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
	        if (cloud.pbc_idx_map[neigh_idx] > idx) {
	          idxs.push_back(idx);
	          neigh_idxs.push_back(cloud.pbc_idx_map[neigh_idx]);
            //std::cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << std::endl;
            //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
            nbs++;
	        }
	      }
        // if the index is not found in the skin map the pair
        // can be stored as is
	      else {
	        idxs.push_back(idx);
	        neigh_idxs.push_back(neigh_idx);
          //std::cout << idx << ", " << neigh_idx << std::endl;
          //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
	        nbs++;
	      }
      }
    }
    }
  }
  for (int i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    for (size_t k = 0; k < nbs; k++) {
      xa = frame.pts[idxs[k]].x;
      xb = frame.pts[neigh_idxs[k]].x;
      ya = frame.pts[idxs[k]].y;
      yb = frame.pts[neigh_idxs[k]].y;
      za = frame.pts[idxs[k]].z;
      zb = frame.pts[neigh_idxs[k]].z;
      if ((xa-xb) < (-xlen*0.5)) {
	      xb = xb - xlen;
      }
      if ((xa-xb) > (xlen*0.5)){
	      xb = xb + xlen;
      }
      if ((ya-yb) < (-ylen*0.5)){
	      yb = yb - ylen;
      }
      if ((ya-yb) > (ylen*0.5)){
	      yb = yb + ylen;
      }
      if ((za-zb) < (-zlen*0.5)){
	      zb = zb - zlen;
      }
      if ((za-zb) > (zlen*0.5)){
	      zb = zb + zlen;
      }
      xdist = std::abs(xa-xb);
      ydist = std::abs(ya-yb);
      zdist = std::abs(za-zb);
      old_avg = avg;
      avg = old_avg + (xdist - old_avg)/n;
      diff_sqrd = diff_sqrd + ((xdist - old_avg)*(xdist - avg));
      n++;
      old_avg = avg;
      avg = old_avg + (ydist - avg)/n;
      diff_sqrd = diff_sqrd + ((ydist - old_avg)*(ydist - avg));
      n++;
      old_avg = avg;
      avg = old_avg + (zdist - avg)/n;
      diff_sqrd = diff_sqrd + ((zdist - old_avg)*(zdist - avg));
      n++;
    }
  }
  std::cout << "average= " << avg << std::endl;
  variance = diff_sqrd/(n-1);
  std::cout << "n pairs= " << n << std::endl;
  std::cout << "variance01= " << variance << std::endl;
  std::cout << "std01= " << pow(variance, 0.5) << std::endl;
}
template <typename num_t>
void variance01kd_r(std::string &filename, simFrame<num_t> &avg_frame, const size_t N,
    const int num_frames, const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 2.5;
  size_t header = 5;
  int n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist_2, ydist_2, zdist_2,
      dist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame;
  //std::cout << "populating point cloud" << std::endl;
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  //std::cout << "building point cloud" << std::endl;
  index.buildIndex();
  size_t num_results = num_nbs;
  int idx;
  int neigh_idx;
  int nbs = 0;
  std::vector<double> idxs;
  std::vector<double> neigh_idxs;
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<num_t> out_dist(num_results);
    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0],
        &out_dist[0]);
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    // downselect neighbors and restore indexes from skin map
    // count number of actual neighbors
    for (size_t k = 0; k < num_results; k++) {
      neigh_idx = ret_index[k];
      if (neigh_idx > idx) {
        // is the neighbor index found in the skin map
	      if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
	        if (cloud.pbc_idx_map[neigh_idx] > idx) {
	          idxs.push_back(idx);
	          neigh_idxs.push_back(cloud.pbc_idx_map[neigh_idx]);
            //std::cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << std::endl;
            //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
            nbs++;
	        }
	      }
        // if the index is not found in the skin map the pair
        // can be stored as is
	      else {
	        idxs.push_back(idx);
	        neigh_idxs.push_back(neigh_idx);
          //std::cout << idx << ", " << neigh_idx << std::endl;
          //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
	        nbs++;
	      }
      }
    }
    }
  }
  for (int i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    for (size_t k = 0; k < nbs; k++) {
      xa = frame.pts[idxs[k]].x;
      xb = frame.pts[neigh_idxs[k]].x;
      ya = frame.pts[idxs[k]].y;
      yb = frame.pts[neigh_idxs[k]].y;
      za = frame.pts[idxs[k]].z;
      zb = frame.pts[neigh_idxs[k]].z;
      if ((xa-xb) < (-xlen*0.5)) {
	      xb = xb - xlen;
      }
      if ((xa-xb) > (xlen*0.5)){
	      xb = xb + xlen;
      }
      if ((ya-yb) < (-ylen*0.5)){
	      yb = yb - ylen;
      }
      if ((ya-yb) > (ylen*0.5)){
	      yb = yb + ylen;
      }
      if ((za-zb) < (-zlen*0.5)){
	      zb = zb - zlen;
      }
      if ((za-zb) > (zlen*0.5)){
	      zb = zb + zlen;
      }
      xdist_2 = pow(std::abs(xa-xb),2.0);
      ydist_2 = pow(std::abs(ya-yb),2.0);
      zdist_2 = pow(std::abs(za-zb),2.0);
      dist = pow(xdist_2 + ydist_2 + zdist_2, 0.5);
      old_avg = avg;
      avg = old_avg + (dist - old_avg)/n;
      diff_sqrd = diff_sqrd + ((dist - old_avg)*(dist - avg));
      n++;
    }
  }
  std::cout << "average= " << avg << std::endl;
  variance = diff_sqrd/(n-1);
  std::cout << "n pairs= " << n << std::endl;
  std::cout << "variance01= " << variance << std::endl;
  std::cout << "std01= " << pow(variance, 0.5) << std::endl;
}
int main() {
  int num_nbs = 8;
  int num_atoms = 16000;
  int num_frames = 5000;
//  std::string filename = "./run_06.dump";
//  std::string filename = "../../../lammps_EF/1K_500.trj";
//  std::string filename = "../../../lammps_EF/100K_5000.trj";
  std::string filename = "../../../lammps_EF/source/testing/1K_500.trj";
//  std::string filename = "./traj_54_0.00002.trj";
//  std::string filename = "../../../lammps_EF/source/testing/run_04.dump";
  resultSet<double> results;
  //NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
  variance00WK<double>(filename, num_atoms, num_frames,results);
  variance01kd_r<double>(filename, results.avg, num_atoms, num_frames, num_nbs);
  std::cout << "variance00= " << results.variance << std::endl;
  std::cout << "std00= " << pow(results.variance,0.5) << std::endl;
  return 0;
}
