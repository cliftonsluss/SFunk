/***********************************************************************
 * Software License Agreement (BSD License)
 *
 * Copyright 2011-2016 Jose Luis Blanco (joseluisblancoc@gmail.com).
 *   All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *************************************************************************/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include "nanoflann.h"


#ifndef TEST_H
#define TEST_H

using namespace nanoflann;

template <class T, class DataSource, typename _DistanceType = T>
struct NN_Adaptor {
  typedef T ElementType;
  typedef _DistanceType DistanceType;

  const DataSource &data_source;

  NN_Adaptor(const DataSource &_data_source)
      : data_source(_data_source) {}

  inline DistanceType evalMetric(const T *a, const size_t b_idx,
                                 size_t size) const {
    DistanceType result = DistanceType();
    for (size_t i = 0; i < size; ++i) {
      const DistanceType diff = a[i] - data_source.kdtree_get_pt(b_idx, i);
      result += ((diff == 0) ? 10000 : diff * diff);
    }
    return result;
  }

  template <typename U, typename V>
  inline DistanceType accum_dist(const U a, const V b, const size_t) const {
    return ((a - b) * (a - b) == 0.0 ? 10000 : (a - b) * (a - b));
  }
};

template <typename T>
struct PointCloud
{
	struct Point
	{
	  T   x,y,z;
	};

	std::vector<Point>  pts;

	struct Atom {
	  int atom_num, atom_type;
	};

	std::vector<Atom>  atms;
  std::vector<int> sphere_list;
  int count;

	struct Bounds {
	  double min,max;
	};

	Bounds Xbounds, Ybounds, Zbounds;

	std::map<int,int> pbc_idx_map;

  // everything below is the stock point cloud from nanoflann
	// Must return the number of data points
	inline size_t kdtree_get_point_count() const { return pts.size(); }

	// Returns the dim'th component of the idx'th point in the class:
	// Since this is inlined and the "dim" argument is typically an immediate value, the
	//  "if/else's" are actually solved at compile time.
	inline T kdtree_get_pt(const size_t idx, const size_t dim) const
	{
		if (dim == 0) return pts[idx].x;
		else if (dim == 1) return pts[idx].y;
		else return pts[idx].z;
	}
	// Optional bounding-box computation: return false to default to a standard bbox computation loop.
	//   Return true if the BBOX was already computed by the class and returned in "bb" so it can be avoided to redo it again.
	//   Look at bb.size() to find out the expected dimensionality (e.g. 2 or 3 for point clouds)
	template <class BBOX>
	bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};


template <typename T>
struct simFrame
{
	struct Point
	{
	  T  x,y,z;
	};
	std::vector<Point>  pts;
	struct Box {
	  T min, max;
	};
	Box xbox, ybox, zbox;
	struct Atom {
	  int atom_num, atom_type;
	};
	std::vector<Atom> atms;
	std::vector<double> flat;
	int num_atoms;
  int index;
};

template <typename T>
struct resultSet
{
  simFrame<T> avg;
  T variance;
};


template <typename T>
// scaled means that fractional coordinates have been converted to
// full value cartesian coordinates
void populatePointCloudPBC(PointCloud<T> &cloud, simFrame<T> &frame,
    double skin, bool scaled = true, const T max_range = 10) {
  // calculate our box widths
  double xlen = frame.xbox.max - frame.xbox.min;
  double ylen = frame.ybox.max - frame.ybox.min;
  double zlen = frame.zbox.max - frame.zbox.min;
  // every atom found outside of these thresholds will be added as skin to
  // opposite side of PointCloud
  double x_thresh_min = frame.xbox.min + skin;
  double x_thresh_max = frame.xbox.max - skin;
  double y_thresh_min = frame.ybox.min + skin;
  double y_thresh_max = frame.ybox.max - skin;
  double z_thresh_min = frame.zbox.min + skin;
  double z_thresh_max = frame.zbox.max - skin;
  cloud.pts.resize(frame.num_atoms);
  if (scaled) {
    // std::cout << "reading full scaled atom positions into point cloud" << std::endl;
  // getNextFrame rescales these so if feeding frames from getNextFrame into
  // populatePointCloudPBC make sure to set scaled = true
    for (int i = 0; i < frame.num_atoms; i++){
      cloud.pts[i].x = frame.pts[i].x;
      cloud.pts[i].y = frame.pts[i].y;
      cloud.pts[i].z = frame.pts[i].z;
    }
  }
  else {
  // this frame was created with fractional coordinates
    for (int i = 0; i < frame.num_atoms; i++) {
      cloud.pts[i].x = (frame.pts[i].x * xlen) + frame.xbox.min;
      cloud.pts[i].y = (frame.pts[i].y * ylen) + frame.ybox.min;
      cloud.pts[i].z = (frame.pts[i].z * zlen) + frame.zbox.min;
    }
  }

  //}
  cloud.count = frame.num_atoms;
  // pad outside of cloud to mimic periodic boundary conditions
  // and create map of added indexes back to original indexes
  // pbc_idx is the index of an atom at it's position in the skin
  // a map of key:value pairs is created with pbc_idx as the key and the
  // original index as the value. Any time an atom is added to the skin more
  // than once it's new pbc_idx is added to pbc_idx_map as the key and the value
  // stored at the last previous pbc_idx is stored as the value.
  int pbc_idx = cloud.count;
  // first we add padding to boundaries in x
  for (int i = 0; i < cloud.count; i++) {
    if (cloud.pts[i].x > x_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x - xlen, cloud.pts[i].y, cloud.pts[i].z});
      cloud.pbc_idx_map[pbc_idx] = i;
      pbc_idx++;
    }
    if (cloud.pts[i].x < x_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x + xlen, cloud.pts[i].y, cloud.pts[i].z});
      cloud.pbc_idx_map[pbc_idx] = i;
      pbc_idx++;
    }
  }
  // we have to reset cloud.count to include the atoms we've added
  cloud.count = pbc_idx;

  // next we pad the y boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].y > y_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x, cloud.pts[i].y - ylen, cloud.pts[i].z});
  // if this point has already been remapped, just map the new reference to
  // the previously mapped point
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	      pbc_idx++;
      }
  // if this point has not already been remapped, then create a new map entry
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
      }
    }
    if (cloud.pts[i].y < y_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x, cloud.pts[i].y + ylen, cloud.pts[i].z});
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	      pbc_idx++;
      }
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
      }
    }
  }

  cloud.count = pbc_idx;
  // finally pad the z boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].z > z_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z - zlen});
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	pbc_idx++;
      }
      else {
	cloud.pbc_idx_map[pbc_idx] = i;
	pbc_idx++;
      }
    }
    if (cloud.pts[i].z < z_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z + zlen});
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	pbc_idx++;
      }
      else {
	cloud.pbc_idx_map[pbc_idx] = i;
	pbc_idx++;
      }
    }
  }
  // std::cout << cloud.pbc_idx_map.size() << " atoms in padding" << std::endl;
  // std::cout << "pbc_idx = " << pbc_idx << std::endl;
}

class RunningStat
{
  public:
    RunningStat() : m_n(0) {}
    void Clear()
    {
      m_n = 0;
    }
    void Push(double x)
    {
      m_n++;
      // See Knuth TAOCP vol 2, 3rd edition, page 232
      if (m_n == 1)
      {
        m_oldM = m_newM = x;
        m_oldS = 0.0;
      }
      else
      {
        m_newM = m_oldM + (x - m_oldM)/m_n;
        m_newS = m_oldS + (x - m_oldM)*(x - m_newM);
        // set up for next iteration
        m_oldM = m_newM;
        m_oldS = m_newS;
      }
    }
    size_t NumDataValues() const
    {
      return m_n;
    }
    double Mean() const
    {
      return (m_n > 0) ? m_newM : 0.0;
    }
    double Variance() const
    {
      return ( (m_n > 1) ? m_newS/(m_n - 1) : 0.0 );
    }
    double StandardDeviation() const
    {
      return sqrt( Variance() );
    }

  private:
    size_t m_n;
    double m_oldM, m_newM, m_oldS, m_newS;
};



class Trajectory;
#include "traj_reader.h"
//class Trajectory {
//  public:
//    Trajectory(std::string &filename, const size_t num_atoms, const size_t header);
//    void getNextFrame(simFrame<double> &frame);//

//  private:
//    std::string filename;
//    size_t num_atoms;
//    size_t header;
//    std::ifstream inputfile;
//};

//Welford method expanded to 3d point cloud
// https://www.johndcook.com/blog/standard_deviation/
template <typename T>
void variance00WK(std::string &filename, int num_atoms, int num_frames,
    int num_skipframes, resultSet<T> &result) {
  double xa, ya, za, xb, yb, zb, xlen, ylen, zlen, old_avgx, old_avgy,
    old_avgz, diff_sqrd, nsamples, variance;
  simFrame<double> frame0;
  simFrame<double> frame;
  //std::vector<double> avg_x, avg_y, avg_z;
  Trajectory traj(filename, num_atoms, 5);
  traj.skipFrames(num_skipframes);
  traj.getNextFrame(frame);
  frame0 = frame;
  result.avg.pts.resize(num_atoms);
  result.avg.num_atoms = num_atoms;

  // for all points in first frame place them as initial values of the average
  // frame.
  // Points are already scaled by Trajectory::getNextFrame method
  for (int j = 0; j < num_atoms; j++) {
    //cout << frame0.pts[j].x << std::endl;
    result.avg.pts[j].x = frame0.pts[j].x;
    result.avg.pts[j].y = frame0.pts[j].y;
    result.avg.pts[j].z = frame0.pts[j].z;
  }
  diff_sqrd = 0;
  double div_3 = 1.0/3.0;

  for (int i = 1; i < num_frames; i++) {
    traj.getNextFrame(frame);
    xlen = frame.xbox.max - frame.xbox.min;
    ylen = frame.ybox.max - frame.ybox.min;
    zlen = frame.zbox.max - frame.zbox.min;
    int k = 0;
    for (int j = 0; j < num_atoms; j++) {
      xa = frame0.pts[j].x;
      xb = frame.pts[j].x;
      ya = frame0.pts[j].y;
      yb = frame.pts[j].y;
      za = frame0.pts[j].z;
      zb = frame.pts[j].z;
      if ((xa - xb) < (-xlen*0.5)){
	      xb = xb - xlen;
      }
      if ((xa - xb) > (xlen*0.5)){
	      xb = xb + xlen;
      }
      if ((ya - yb) < (-ylen*0.5)){
	      yb = yb - ylen;
      }
      if ((ya - yb) > (ylen*0.5)){
	      yb = yb + ylen;
      }
      if ((za - zb) < (-zlen*0.5)){
	      zb = zb - zlen;
      }
      if ((za - zb) > (zlen*0.5)){
	      zb = zb + zlen;
      }
      old_avgx = result.avg.pts[j].x;
      old_avgy = result.avg.pts[j].y;
      old_avgz = result.avg.pts[j].z;
      result.avg.pts[j].x = old_avgx + (xb - result.avg.pts[j].x)/i;
      result.avg.pts[j].y = old_avgy + (yb - result.avg.pts[j].y)/i;
      result.avg.pts[j].z = old_avgz + (zb - result.avg.pts[j].z)/i;
      diff_sqrd = diff_sqrd + ((xb - old_avgx)*(xb - result.avg.pts[j].x)
	        + (yb - old_avgy)*(yb - result.avg.pts[j].y)
          + (zb - old_avgz)*(zb - result.avg.pts[j].z));

      if (result.avg.pts[j].x < result.avg.xbox.min) {
	result.avg.xbox.min = result.avg.pts[j].x;
      }
      if (result.avg.pts[j].x > result.avg.xbox.max) {
	result.avg.xbox.max = result.avg.pts[j].x;
      }
      if (result.avg.pts[j].y < result.avg.ybox.min) {
	result.avg.ybox.min = result.avg.pts[j].y;
      }
      if (result.avg.pts[j].y > result.avg.ybox.max) {
	result.avg.ybox.max = result.avg.pts[j].y;
      }
      if (result.avg.pts[j].z < result.avg.zbox.min) {
	result.avg.zbox.min = result.avg.pts[j].z;
      }
      if (result.avg.pts[j].z > result.avg.zbox.max) {
	result.avg.zbox.max = result.avg.pts[j].z;
      }
    }
  }
  // what is this bit here? ah, we're returning to fractional coordinates
  // double xfac = 1.0/(result.avg.xbox.max - result.avg.xbox.min);
  // double yfac = 1.0/(result.avg.ybox.max - result.avg.ybox.min);
  // double zfac = 1.0/(result.avg.zbox.max - result.avg.zbox.min);
  // for (int i = 0; i < num_atoms; i++) {
  //   result.avg.pts[i].x = (result.avg.pts[i].x - result.avg.xbox.min)*xfac;
  //   result.avg.pts[i].y = (result.avg.pts[i].y - result.avg.ybox.min)*yfac;
  //   result.avg.pts[i].z = (result.avg.pts[i].z - result.avg.zbox.min)*zfac;
  // }

  diff_sqrd = diff_sqrd/num_atoms;
  nsamples = num_frames*3.0;
  result.variance = diff_sqrd/(nsamples-1);


  return;
}

// below copied from original traj.cpp

// function to calculate mean and standard deviation given a set of data in the
// std::vector "data" returns a std::vector {mean, standard deviation}
// std::vector<double> vectorStats(std::vector<double> data) {
//   double sum = std::accumulate(data.begin(), data.end(), 0.0);
//   double mean = sum/data.size();
//   std::vector<double> diff(data.size());
//   transform(data.begin(), data.end(), diff.begin(),
//       bind2nd(std::minus<double>(), mean));
//   double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
//   double stdev = sqrt(sq_sum / data.size());
//   std::vector<double> out = {mean, stdev};
//   return out;
// }
template <typename num_t>
void NNeighbors(std::string &filename, const size_t N, const int num_frames,
    const int num_nbs) {
  PointCloud<num_t> cloud;
  double skin = 4.0;
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
  double skin = 6.0;
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
  // grab our neighbors and apply minimum image convention
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
    const int num_frames, const int num_skipframes, const int num_nbs,
    size_t &nbs_found, double &variance01, std::string &outfile,
    int dump=0, float error=0.0001) {
  PointCloud<num_t> cloud;
  if (!outfile.empty()) {
    std::ofstream var_out;
    var_out.open(outfile);
    var_out << "npairs\t   variance\n";
    var_out.close();
  }

  double skin = 4.0;
  size_t header = 5;
  size_t n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist_2, ydist_2, zdist_2,
      dist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  double var_check;
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
  size_t min_count = 100;
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
  // std::cout << "neighbor count " << nbs << std::endl;
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

      var_check = rs.Variance();
      rs.Push(dist);
      // if (rs.NumDataValues() > min_count){
      //   if (std::abs(var_check - rs.Variance()) < error){
      //     nbs_found = rs.NumDataValues();
      //     variance01 = rs.Variance();
      //     std::cout << "Achieved minimum error with " << nbs_found << " pairs\n";
      //     std::cout << "variance = " << rs.Variance() << "\n";
      //     return;
      //   }
      // }


      if (dump > 0) {
        if (rs.NumDataValues() % dump == 0) {
          std::ofstream var_out;
          var_out.open(outfile, std::ofstream::app);
          var_out << rs.NumDataValues() << " " <<  rs.Variance() << "\n";
          var_out.close();
        }
      }

      // variance_vector.push_back(rs.Variance());

    }
  }
  // std::cout << "average= " << avg << std::endl;
  // std::cout << "average= " << rs.Mean() << std::endl;
  // variance = diff_sqrd/(n-1);
  // std::cout << "n pairs= " << n << std::endl;
  nbs_found = rs.NumDataValues();
  // std::cout << "number of pairs= " << nbs_found << std::endl;
  // std::cout << "variance01= " << variance << std::endl;
  variance01 = rs.Variance();
  // std::cout << "variance01= " << variance01 << std::endl;
  // std::cout << "std01= " << pow(variance, 0.5) << std::endl;


}




#endif
