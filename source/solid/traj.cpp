#include "nanoflann.h"
#include "read_config.h"
#include "traj_reader.h"
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
  frame.index = 0;
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
    const int num_frames, const int num_skipframes, const int num_nbs) {
// create PointCloud object
  PointCloud<num_t> cloud;
  // skin defines amount of padding to add to outsides of simulation cube
  // units are Angstroms
  // Instead of the minimum image criterium which won't work for KDtree
  // we replicate enough of the region from PBC necessary to perform calculations
  double skin = 5.0;
  size_t header = 5;
  long n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist, zdist,
      variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame;
  //std::cout << "populating point cloud" << std::endl;
  // is avg_frame scaled or no?
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  // dimensionality = 3, DatasetAdaptor = cloud, max leaf count = 10
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  //std::cout << "building point cloud" << std::endl;
  index.buildIndex();
  size_t num_results = num_nbs;
  long idx;
  long neigh_idx;
  long nbs = 0;
  std::vector<double> idxs;
  std::vector<double> neigh_idxs;
  // for all atoms
  for (int j = 0; j < N; j++)
  {
    idx = j;
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<num_t> out_dist(num_results);
  // pass output vectors by reference
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
  // armed with a neighbor list based on an average frame, we calculate the
  // neighbor distances
  RunningStat rs;
  traj.skipFrames(num_skipframes);
  for (int i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
  // grab our neighbors and apply minimum image criterium
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
  // finally get the absolute value of the difference
      xdist = std::abs(xa-xb);
      ydist = std::abs(ya-yb);
      zdist = std::abs(za-zb);
  // below is an implementation of the Welford method of running statistics
  // extended to 3 dimensions, all dimensions are summed into one dimension
      // old_avg = avg;
      // avg = old_avg + (xdist - old_avg)/n;
      // diff_sqrd = diff_sqrd + ((xdist - old_avg)*(xdist - avg));
      // n++;
      // old_avg = avg;
      // avg = old_avg + (ydist - old_avg)/n;
      // diff_sqrd = diff_sqrd + ((ydist - old_avg)*(ydist - avg));
      // n++;
      // old_avg = avg;
      // avg = old_avg + (zdist - old_avg)/n;
      // diff_sqrd = diff_sqrd + ((zdist - old_avg)*(zdist - avg));
      // n++;
      rs.Push(xdist);
      rs.Push(ydist);
      rs.Push(zdist);

    }
  }
  // std::cout << "average= " << avg << std::endl;
  std::cout << "average= " << rs.Mean() << std::endl;
  // variance = diff_sqrd/(n-1);
  // std::cout << "n pairs= " << n << std::endl;
  // std::cout << "variance01= " << variance << std::endl;
  std::cout << "variance01= " << rs.Variance() << std::endl;
  std::cout << "std01= " << pow(rs.Variance(), 0.5) << std::endl;
}
template <typename num_t>
void variance01kd_r(std::string &filename, simFrame<num_t> &avg_frame, const size_t N,
    const int num_frames, const int num_skipframes, const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 4.0;
  size_t header = 5;
  size_t n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist_2, ydist_2, zdist_2,
      dist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame;
  //std::cout << "populating point cloud" << std::endl;
  // we are using the average positions of the atoms from all of the frames to
  // generate a nearest neighbor list
  // populate point cloud with average frame
  populatePointCloudPBC(cloud, avg_frame, skin);
  // set up KDtree adaptor
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  // create an index for our adaptor
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
  //std::cout << "building point cloud" << std::endl;
  index.buildIndex();
  size_t num_results = num_nbs;
  size_t idx;
  size_t neigh_idx;
  size_t nbs = 0;
  std::vector<double> idxs;
  std::vector<double> neigh_idxs;
  // we will use each atom as a query point to find nearest neighbors
  for (size_t j = 0; j < N; j++)
  {
    idx = j;
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z };
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<num_t> out_dist(num_results);
    // knnSearch returns size_t resultSet.size() num_results is the number of
    // neighbors we requested for this U.C. type
    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0],
        &out_dist[0]);
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    // downselect neighbors and restore indexes from skin map
    // count number of actual neighbors
    // std::cout << "atom index, number of neighbors= " << j << ", " << num_results << std::endl;
    for (size_t k = 0; k < num_results; k++) {
      // std::cout << "neighbor count " << nbs << std::endl;
      neigh_idx = ret_index[k];
      // for any given neighbor index (neigh_idx) that is greater than the
      // current atom's index (idx) we want to add it to the list of neighbor
      // indexs (neigh_idxs) the greater than condition ensures we don't include
      // duplicate pairs
      if (neigh_idx > idx) {
        // is the neighbor index found in the skin map
	      if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) {
	        if (cloud.pbc_idx_map[neigh_idx] > idx) {
	          idxs.push_back(idx);
            // std::cout << "index " << idx << std::endl;
	          neigh_idxs.push_back(cloud.pbc_idx_map[neigh_idx]);
            // std::cout << "neighbor index " << neigh_idx << std::endl;
            //std::cout << idx << ", " << cloud.pbc_idx_map[neigh_idx] << std::endl;
            //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
            nbs++;
	        }
	      }
        // if the index is not found in the skin map the pair
        // can be stored as is
	      else {
	        idxs.push_back(idx);
          // std::cout << "index " << idx << std::endl;
	        neigh_idxs.push_back(neigh_idx);
          // std::cout << "neighbor index " << neigh_idx << std::endl;
          //std::cout << idx << ", " << neigh_idx << std::endl;
          //std::cout << idxs[nbs] << ", " << neigh_idxs[nbs] << std::endl;
	        nbs++;
	      }
      }
    }
    }
  }
  RunningStat rs;
  std::cout << "neighbor count " << nbs << std::endl;
  traj.skipFrames(num_skipframes);
  for (size_t i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    // std::cout << "frame number = " << i << std::endl;
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    // std::cout << "nbs = " << nbs << std::endl;
    for (size_t k = 0; k < nbs; k++) {
      // grab our reference point x,y,z coordinates and coordinates of it's
      // neighbors
      xa = frame.pts[idxs[k]].x;
      xb = frame.pts[neigh_idxs[k]].x;
      ya = frame.pts[idxs[k]].y;
      yb = frame.pts[neigh_idxs[k]].y;
      za = frame.pts[idxs[k]].z;
      zb = frame.pts[neigh_idxs[k]].z;
      // apply minimum image criteria and calculate neighbor distances in xyz
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
      // we will calculate the average sum of the squared component distances
      // and then take the square root of the final result
      xdist_2 = pow((xa-xb),2.0);
      ydist_2 = pow((ya-yb),2.0);
      zdist_2 = pow((za-zb),2.0);
      dist = pow(xdist_2 + ydist_2 + zdist_2, 0.5);
      // std::cout << "atom " <<

      // rs.Push(dist);

      old_avg = avg;
      avg = old_avg + (dist - old_avg)/n;
      // diff_sqrd/(n-1) will be our variance
      diff_sqrd = diff_sqrd + (dist - old_avg)*(dist - avg);
      n++;
    }
  }
  std::cout << "average= " << avg << std::endl;
  // std::cout << "average= " << rs.Mean() << std::endl;
  variance = diff_sqrd/(n-1);
  std::cout << "n pairs= " << n << std::endl;
  std::cout << "variance01= " << variance << std::endl;
  // std::cout << "variance01= " << rs.Variance() << std::endl;
  // std::cout << "std01= " << pow(variance, 0.5) << std::endl;
}

int main(int argc, char *argv[]) {
  // read parameters from json file
  std::map<std::string,std::string> config_map;
  std::string config_file = argv[1];
  std::cout << "Reading configuration file..." << std::endl;
  Read_config config(config_file, config_map);
  config.get_config();
  int num_nbs = config.j["neighbors"];
  std::cout << num_nbs << " nearest neighbors" << std::endl;
  int num_atoms = config.j["atoms"];
  std::cout << num_atoms << " atoms" << std::endl;
  int num_frames = config.j["frames"];
  std::cout << num_frames << " frames" << std::endl;
  int num_skipframes = config.j["skipframes"];
  std::cout << "Skipping " << num_skipframes << " frames" << std::endl;
  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;


//lines below should all remove // to uncomment
// create a resultSet obect to house our results
resultSet<double> results;
//NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
// only results.avg is needed from variance00WK at thi point, that functionality
// will be broken out into a more compact method at a later date
variance00WK<double>(datafile, num_atoms, num_frames, num_skipframes, results);
// values below not needed anymore but still being reported for reference
std::cout << "variance00= " << results.variance << std::endl;
std::cout << "std00= " << pow(results.variance,0.5) << std::endl;

variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames, num_skipframes, num_nbs);
// variance01kd<double>(datafile, results.avg, num_atoms, num_frames, num_skipframes, num_nbs);

return 0;
}
