#include "nanoflann.h"
#include "test.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>
using namespace std;
using namespace nanoflann;
// function to calculate mean and standard deviation given a set of data in the
// vector "data" returns a vector {mean, standard deviation}
vector<double> VectorStats(vector<double> data) {
  double sum = accumulate(data.begin(), data.end(), 0.0);
  double mean = sum/data.size();
  vector<double> diff(data.size());
  transform(data.begin(), data.end(), diff.begin(),
      bind2nd(minus<double>(), mean));
  double sq_sum = inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
  double stdev = sqrt(sq_sum / data.size());
  vector<double> out = {mean, stdev};
  return out;
}
template <typename num_t>
void NNeighbors(string &filename, const size_t N, const int num_frames,
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
    vector<size_t> idxs;
    vector<size_t> neigh_idxs;
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    for (int j = 0; j < N; j++)
    {
      idx = j;
      const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
          cloud.pts[idx].z };
      {
      vector<size_t> ret_index(num_results);
      vector<num_t> out_dist(num_results);
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
              //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << endl;
              nbs++;
  	        }
  	      }
          // if the index is not found in the skin map the pair
          // can be stored as is
  	      else {
  	        idxs.push_back(idx);
  	        neigh_idxs.push_back(neigh_idx);
            //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << endl;
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
        xdist = abs(xa-xb);
        ydist = abs(ya-yb);
        zdist = abs(za-zb);
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
  cout << "average= " << avg << endl;
  variance = diff_sqrd/(n-1);
  cout << "n pairs= " << n << endl;
  cout << "variance01= " << variance << endl;
  cout << "std01= " << pow(variance, 0.5) << endl;
}

template <typename num_t>
void variance01kd(string &filename, simFrame<num_t> &avg_frame, const size_t N,
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
  //cout << "populating point cloud" << endl;
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  //cout << "building point cloud" << endl;
  index.buildIndex();
  size_t num_results = num_nbs;
  int idx;
  int neigh_idx;
  int nbs = 0;
  vector<double> idxs;
  vector<double> neigh_idxs;
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    vector<size_t> ret_index(num_results);
    vector<num_t> out_dist(num_results);
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
            //cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << endl;
            //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << endl;
            nbs++;
	        }
	      }
        // if the index is not found in the skin map the pair
        // can be stored as is
	      else {
	        idxs.push_back(idx);
	        neigh_idxs.push_back(neigh_idx);
          //cout << idx << ", " << neigh_idx << endl;
          //cout << idxs[nbs] << ", " << neigh_idxs[nbs] << endl;
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
      xdist = abs(xa-xb);
      ydist = abs(ya-yb);
      zdist = abs(za-zb);
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
  cout << "average= " << avg << endl;
  variance = diff_sqrd/(n-1);
  cout << "n pairs= " << n << endl;
  cout << "variance01= " << variance << endl;
  cout << "std01= " << pow(variance, 0.5) << endl;
}

int main() {
  int num_nbs = 8;
  int num_atoms = 16000;
  int num_frames = 500;
//  string filename = "./run_06.dump";
  string filename = "../../../lammps_EF/1K_500.trj";
//  string filename = "../../../lammps_EF/100K_5000.trj";
//  string filename = "../../../lammps_EF/source/testing/1K_5000.trj";
//  string filename = "./traj_54_0.00002.trj";
  resultSet<double> results;
  NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
  variance00WK<double>(filename, num_atoms, num_frames,results);
  variance01kd<double>(filename, results.avg, num_atoms, num_frames, num_nbs);
  cout << "variance00= " << results.variance << endl;
  cout << "std00= " << pow(results.variance,0.5) << endl;
  return 0;
}
