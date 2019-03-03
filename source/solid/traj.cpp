#include "nanoflann.h"
#include "test.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <numeric>
#include <algorithm>

using namespace std;
using namespace nanoflann;

// function to calculate mean and standard deviation given a set of data in the vector
// "data" returns a vector {mean, standard deviation}
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
void NNeighbors(string &filename, const size_t N, const int num_frames, const int num_nbs) {

  PointCloud<num_t> cloud;

  //string filename = "./iron_test.dump"; 
  double skin = 2.5;
  size_t header = 5;
  size_t tot_frames = 5;
  size_t frame_num = 3;
  int n = 1;
  double xb, yb, zb, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  

  vector<double> xdist, ydist, zdist, L2dist; 
  
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame; 
  for (int i = 0; i < num_frames; i++) { 
    
    
    traj.getNextFrame(frame);
    populatePointCloudPBC(cloud, frame, skin);
  
    //construct a kd-tree index
    typedef KDTreeSingleIndexAdaptor<
      L2_Simple_Adaptor<num_t, PointCloud<num_t> >,
      //NN_Adaptor<num_t, PointCloud<num_t> >,
      PointCloud<num_t>,
      3
      > my_kd_tree_t;

    my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
    index.buildIndex();

    size_t num_results = num_nbs;
  
    int idx;
    int neigh_idx;

    for (int j = 0; j < N; j++)
    { 
      idx = j; 
      const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y, cloud.pts[idx].z };
    
      { 
      vector<size_t> ret_index(num_results);
      vector<num_t> out_dist(num_results);

      num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist[0]);

      ret_index.resize(num_results);
      out_dist.resize(num_results); 
      for (size_t k = 0; k < num_results; k++) {
        neigh_idx = ret_index[k];
        if (neigh_idx > idx) {
	  xb = abs(cloud.pts[idx].x - cloud.pts[neigh_idx].x);
          yb = abs(cloud.pts[idx].y - cloud.pts[neigh_idx].y);
          zb = abs(cloud.pts[idx].z - cloud.pts[neigh_idx].z);
  
	  old_avg = avg;
	  avg = old_avg + (xb - old_avg)/n;
	  diff_sqrd = diff_sqrd + ((xb - old_avg)*(xb - avg));
	  n++;
	  old_avg = avg;
	  avg = old_avg + (yb - avg)/n;
	  diff_sqrd = diff_sqrd + ((yb - old_avg)*(yb - avg));
	  n++;
	  old_avg = avg;
	  avg = old_avg + (zb - avg)/n;
	  diff_sqrd = diff_sqrd + ((zb - old_avg)*(zb - avg));
	  n++;

          //xdist.push_back(abs(cloud.pts[idx].x - cloud.pts[neigh_idx].x));
          //ydist.push_back(abs(cloud.pts[idx].y - cloud.pts[neigh_idx].y));
          //zdist.push_back(abs(cloud.pts[idx].z - cloud.pts[neigh_idx].z));
	  //L2dist.push_back(out_dist[k]);
        }

      }
      }
    }
  }
  variance = diff_sqrd/(n-1);
  //cout << "n pairs= " << n << endl;
  cout << "variance01= " << variance << endl;
  cout << "std01= " << pow(variance, 0.5) << endl;
  //vector<double> xstats = VectorStats(xdist);
  //cout << "x mean= " << xstats[0] << endl;
  //cout << "x stdev= " << xstats[1] << endl;
  //cout << "x len= " << xdist.size() << endl;
  //vector<double> ystats = VectorStats(ydist);
  //cout << "y mean= " << ystats[0] << endl;
  //cout << "y stdev= " << ystats[1] << endl;
  //vector<double> zstats = VectorStats(zdist);
  //cout << "z mean= " << zstats[0] << endl;
  //cout << "z stdev= " << zstats[1] << endl;
}

template <typename num_t>
void variance01kd(string &filename, simFrame<num_t> &avg_frame, const size_t N, const int num_frames, const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 2.5;
  size_t header = 5;
  int n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist, zdist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<num_t> frame; 
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<num_t, PointCloud<num_t> >,
   PointCloud<num_t>,
   3
   > my_kd_tree_t;
  my_kd_tree_t index(3, cloud, KDTreeSingleIndexAdaptorParams(10) );
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
    const num_t query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y, cloud.pts[idx].z };
    {
    vector<size_t> ret_index(num_results);
    vector<num_t> out_dist(num_results);
    num_results = index.knnSearch(&query_pt[0], num_results, &ret_index[0], &out_dist[0]);
    ret_index.resize(num_results);
    out_dist.resize(num_results);
    //cout << "atom num= " << idx << ":" << endl;
    for (size_t k = 0; k < num_results; k++) {
      neigh_idx = ret_index[k];
      if (neigh_idx > idx) {
	//cout << "old pair = " << idx << ", " << neigh_idx << endl;
	if (cloud.pbc_idx_map.find(neigh_idx) != cloud.pbc_idx_map.end()) { 
	  if (cloud.pbc_idx_map[neigh_idx] > idx) {
	    //cout << "neigh changed" << endl; 
	    idxs.push_back(idx);
	    neigh_idxs.push_back(cloud.pbc_idx_map[neigh_idx]);
	    //cout << "idx= " << idx << "   neigh= " << neigh_idx << endl;
	    nbs++;
	  }
	}
	else {
	  idxs.push_back(idx);
	  neigh_idxs.push_back(neigh_idx);
	  //cout << "idx= " << idx << "   neigh= " << neigh_idx << endl;
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
      //cout << "frame= " << i << ", idx= " << idxs[k] << ", neigh_idx= " << neigh_idxs[k] << endl;
      //cout << "xdist= " << xdist << endl;
      //cout << "ydist= " << ydist << endl;
      //cout << "zdist= " << zdist << endl;
      //cout << "Rdist= " << pow(xdist*xdist + ydist*ydist + zdist*zdist, 0.5) << endl;
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
//  string filename = "../../100K_5000.trj";
  string filename = "../../100K_5000.trj";
  resultSet<double> results;
  variance00WK<double>(filename, num_atoms, num_frames,results);
  NNeighbors<double>(filename, num_atoms, num_frames, num_nbs);
  //variance01kd<double>(filename, results.avg, num_atoms, num_frames, num_nbs);
  cout << "variance00= " << results.variance << endl;
  cout << "std00= " << pow(results.variance,0.5) << endl;
  return 0;
}
