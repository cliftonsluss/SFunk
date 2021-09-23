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
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <sstream>
#include <map>
#include <stdexcept>
#include "nanoflann.h"
#include "neighbor_list_generator.h"
#include "PBCPointCloud_generator.h"
#include "traj_reader.h"
#include "structures.h"

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


// scaled means that fractional coordinates have been converted to
// full value cartesian coordinates
template <typename T>
void cloud2frame(PointCloud<T> &cloud, simFrame<T> &frame) {
  size_t num_atoms = cloud.pts.size();
  std::cout << "num_atoms = " << num_atoms << std::endl;
  double xmin = 0;
  double xmax = 0;
  double ymin = 0;
  double ymax = 0;
  double zmin = 0;
  double zmax = 0;

  // std::cout << "set initial maxes\n";
  // std::cout << "cloud x665 = " << cloud.pts[664].x << std::endl;
  // std::cout << "frame x665 = " << frame.pts[664].x << std::endl;

  for (size_t i = 0; i < num_atoms; i++){
    frame.pts[i].x = cloud.pts[i].x;
    frame.pts[i].y = cloud.pts[i].y;
    frame.pts[i].z = cloud.pts[i].z;

    if (frame.pts[i].x < xmin) {
      xmin = frame.pts[i].x;
    }
    if (frame.pts[i].x > xmax) {
      xmax = frame.pts[i].x;
    }
    if (frame.pts[i].y < ymin) {
      ymin = frame.pts[i].y;
    }
    if (frame.pts[i].y > ymax) {
      ymax = frame.pts[i].y;
    }
    if (frame.pts[i].z < zmin) {
      zmin = frame.pts[i].z;
    }
    if (frame.pts[i].z > zmax) {
      zmax = frame.pts[i].z;
    }
    // std::cout << i << std::endl;
  }
  // std::cout << "frame x665 = " << frame.pts[664].x << std::endl;
  // std::cout << "calculated frame maxes\n";

  frame.box.xmin = xmin - 2;
  frame.box.xmax = xmax + 2;
  frame.box.ymin = ymin - 2;
  frame.box.ymax = ymax + 2;
  frame.box.zmin = zmin - 2;
  frame.box.zmax = zmax + 2;

  frame.box.xlen = xmax - xmin;
  frame.box.ylen = ymax - ymin;
  frame.box.zlen = zmax - zmin;

  frame.num_atoms = num_atoms;
  std::cout << "converted cloud to frame\n";
}

template <typename T>
void populatePointCloudPBC(PointCloud<T> &cloud, simFrame<T> &frame,
    double skin, bool scaled = true, const T max_range = 10) {

  double xlen = frame.box.xlen;
  double ylen = frame.box.ylen;
  double zlen = frame.box.zlen;

  // std::cout << "xlen = " << frame.box.xlen << std::endl;
  // std::cout << "ylen = " << frame.box.ylen << std::endl;
  // std::cout << "zlen = " << frame.box.zlen << std::endl;


  // xlen = frame.box.xmax - frame.box.xmin;
  // ylen = frame.box.ymax - frame.box.ymin;
  // zlen = frame.box.zmax - frame.box.zmin;

  // std::cout << "xlen = " << frame.box.xlen << std::endl;
  // std::cout << "ylen = " << frame.box.ylen << std::endl;
  // std::cout << "zlen = " << frame.box.zlen << std::endl;

  // every atom found outside of these thresholds will be added as skin to
  // opposite side of PointCloud
  double x_thresh_min = frame.box.xmin + skin;
  double x_thresh_max = frame.box.xmax - skin;
  double y_thresh_min = frame.box.ymin + skin;
  double y_thresh_max = frame.box.ymax - skin;
  double z_thresh_min = frame.box.zmin + skin;
  double z_thresh_max = frame.box.zmax - skin;
  cloud.pts.resize(frame.num_atoms);
  if (scaled) {
    // std::cout << "reading full scaled atom positions into point cloud" << std::endl;
  // getNextFrame rescales these so if feeding frames from getNextFrame into
  // populatePointCloudPBC make sure to set scaled = true
    for (int i = 0; i < frame.num_atoms; i++){
      cloud.pts[i].x = frame.pts[i].x;
      cloud.pts[i].y = frame.pts[i].y;
      cloud.pts[i].z = frame.pts[i].z;
      cloud.pbc_idx_map[i] = i;
    }
  }
  else {
  // this frame was created with fractional coordinates
    for (int i = 0; i < frame.num_atoms; i++) {
      cloud.pts[i].x = (frame.pts[i].x * frame.box.xlen) + frame.box.xmin;
      cloud.pts[i].y = (frame.pts[i].y * frame.box.ylen) + frame.box.ymin;
      cloud.pts[i].z = (frame.pts[i].z * frame.box.zlen) + frame.box.zmin;
      cloud.pbc_idx_map[i] = i;
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
  std::cout << "pbc_idx start = " << pbc_idx << std::endl;
  int diff_idx = pbc_idx;
  // first we add padding to boundaries in x
  for (int i = 0; i < cloud.count; i++) {
    if (cloud.pts[i].x > x_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x - xlen,
                           cloud.pts[i].y,
                           cloud.pts[i].z});
      cloud.pbc_idx_map[pbc_idx] = i;

      pbc_idx++;
      // std::cout << "x of index " << i << " modified by subtracting " << xlen
      // << " and added as new index " << pbc_idx-1 << std::endl;
    }
    if (cloud.pts[i].x < x_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x + xlen,
                           cloud.pts[i].y,
                           cloud.pts[i].z});
      cloud.pbc_idx_map[pbc_idx] = i;
      if (i==93243){
        std::cout << "point 93243 mapped to " << pbc_idx << " at " << "x = "
        << cloud.pts[pbc_idx].x<< ", y = " <<
        cloud.pts[pbc_idx].y << ", z = " << cloud.pts[pbc_idx].z << std::endl;
      }
      pbc_idx++;
      // std::cout << "x of index " << i << " modified by adding " << xlen
      // << " and added as new index " << pbc_idx-1 << std::endl;
    }
  }
  std::cout << pbc_idx-diff_idx << " atoms added in x dir\n";
  diff_idx = pbc_idx;
  // we have to reset cloud.count to include the atoms we've added
  cloud.count = pbc_idx;

  // next we pad the y boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].y > y_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y - ylen,
                           cloud.pts[i].z});
     if (i==93243){
       std::cout << "point 93243 y_mapped to " << pbc_idx << " at " << "x = "
       << cloud.pts[pbc_idx].x << ", y = " <<
       cloud.pts[pbc_idx].y << ", z = " << cloud.pts[pbc_idx].z << std::endl;
     }
  // if this point has already been remapped, ensure the new reference is
  // mapped to an original index
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	      pbc_idx++;
        // std::cout << "y of index " << i << " modified by subtracting " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }

  // if this point has not already been remapped, then create a new map entry
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
        // std::cout << "y of index " << i << " modified by subtracting " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
    }
    if (cloud.pts[i].y < y_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y + ylen,
                           cloud.pts[i].z});
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	      pbc_idx++;
        // std::cout << "y of index " << i << " modified by adding " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
        // std::cout << "y of index " << i << " modified by adding " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
    }
  }
  std::cout << pbc_idx-diff_idx << " atoms added in y dir\n";
  diff_idx = pbc_idx;

  cloud.count = pbc_idx;
  // finally pad the z boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].z > z_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y,
                           cloud.pts[i].z - zlen});

       if (i==116571){
         std::cout << "point 116571 z_mapped to " << pbc_idx << " at " << "x = "
         << cloud.pts[pbc_idx].x << ", y = " <<
         cloud.pts[pbc_idx].y << ", z = " << cloud.pts[pbc_idx].z << std::endl;
       }
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
        if (pbc_idx==132663){
          std::cout << "132663 mapped to " << cloud.pbc_idx_map[i] << std::endl;
        }
	      pbc_idx++;
        // std::cout << "z of index " << i << " modified by subtracting " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
        // std::cout << "z of index " << i << " modified by subtracting " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
    }
    if (cloud.pts[i].z < z_thresh_min) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y,
                           cloud.pts[i].z + zlen});
      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
	      pbc_idx++;
        // std::cout << "z of index " << i << " modified by adding " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
      else {
	      cloud.pbc_idx_map[pbc_idx] = i;
	      pbc_idx++;
        // std::cout << "z of index " << i << " modified by adding " << xlen
        // << " and added as new index " << pbc_idx-1 << std::endl;
      }
    }
  }
  std::cout << pbc_idx-diff_idx << " atoms added in z dir\n";
  std::cout << cloud.pbc_idx_map.size() << " atoms in padding" << std::endl;
  std::cout << "pbc_idx = " << pbc_idx << std::endl;

}

template <typename T>
void populatePointCloudSkin(PointCloud<T> &cloud, simFrame<T> &frame,
    double skin, bool scaled = true, const T max_range = 10) {

  double xlen = frame.box.xlen;
  double ylen = frame.box.ylen;
  double zlen = frame.box.zlen;

  // every atom found outside of these thresholds will be added as skin to
  // opposite side of PointCloud
  double x_thresh_min = frame.box.xmin + skin;
  double x_thresh_max = frame.box.xmax - skin;
  double y_thresh_min = frame.box.ymin + skin;
  double y_thresh_max = frame.box.ymax - skin;
  double z_thresh_min = frame.box.zmin + skin;
  double z_thresh_max = frame.box.zmax - skin;

  cloud.pts.resize(frame.num_atoms);
  if (scaled) {
    // std::cout << "reading full scaled atom positions into point cloud" << std::endl;
  // getNextFrame rescales these so if feeding frames from getNextFrame into
  // populatePointCloudPBC make sure to set scaled = true
    for (int i = 0; i < frame.num_atoms; i++){
      cloud.pts[i].x = frame.pts[i].x;
      cloud.pts[i].y = frame.pts[i].y;
      cloud.pts[i].z = frame.pts[i].z;
      cloud.pbc_idx_map[i] = i;
    }
  }
  else {
  // this frame was created with fractional coordinates
    for (int i = 0; i < frame.num_atoms; i++) {
      cloud.pts[i].x = (frame.pts[i].x * frame.box.xlen) + frame.box.xmin;
      cloud.pts[i].y = (frame.pts[i].y * frame.box.ylen) + frame.box.ymin;
      cloud.pts[i].z = (frame.pts[i].z * frame.box.zlen) + frame.box.zmin;
      cloud.pbc_idx_map[i] = i;
    }
  }

  //}
  cloud.count = frame.num_atoms;
  // cloud.list.push_back(27);
  // for (int i = 0; i < cloud.list.size(); i++){
  //   std::cout << "neighbors to be found for the following indexes\n"
  //   << cloud.list[i] << " " << cloud.pts[cloud.list[i]].x << ", " <<
  //   cloud.pts[cloud.list[i]].y << ", " <<
  //   cloud.pts[cloud.list[i]].z << "\n";
  // }
  for (int i = 0; i < cloud.count; i++) {
    // cloud.pts.push_back({cloud.pts[i].x,
    //                      cloud.pts[i].y,
    //                      cloud.pts[i].z});
    if (cloud.pts[i].x < x_thresh_max && cloud.pts[i].x > x_thresh_min) {
      if (cloud.pts[i].y < y_thresh_max && cloud.pts[i].y > y_thresh_min) {
        if (cloud.pts[i].z < z_thresh_max && cloud.pts[i].z > z_thresh_min) {
          // std::cout << i << ", " << cloud.pts[i].x << ", " << cloud.pts[i].y <<
          // ", " << cloud.pts[i].z << std::endl;
          cloud.list.push_back(i);
        }
      }
    }
  }
};

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

//Welford method expanded to 3d point cloud
// https://www.johndcook.com/blog/standard_deviation/
template <typename T>
void variance00WK(std::string &filename, int num_atoms, int num_frames,
    int num_skipframes, resultSet<T> &result) {
  double xa, ya, za, xb, yb, zb, xlen, ylen, zlen, old_avgx, old_avgy,
    old_avgz, diff_sqrd, nsamples, variance;
  simFrame<double> frame0;
  simFrame<double> frame;
  Trajectory traj(filename, num_atoms, 5);
  traj.skipFrames(num_skipframes);
  traj.getNextFrame(frame);
  frame0 = frame;
  result.avg.pts.resize(num_atoms);
  result.avg.num_atoms = num_atoms;
  result.avg.box.xmin = frame.box.xmin;
  result.avg.box.xmax = frame.box.xmax;
  result.avg.box.ymin = frame.box.ymin;
  result.avg.box.ymax = frame.box.ymax;
  result.avg.box.zmin = frame.box.zmin;
  result.avg.box.zmax = frame.box.zmax;
  xlen = frame.box.xlen;
  ylen = frame.box.ylen;
  zlen = frame.box.zlen;


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

  RunningStat rs;
  for (int i = 1; i < num_frames; i++) {
    traj.getNextFrame(frame);

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
    }
  }

  for (int j = 0; j < num_atoms; j++){
    if (result.avg.pts[j].x < result.avg.box.xmin) {
result.avg.pts[j].x = result.avg.pts[j].x + xlen;
    }
    if (result.avg.pts[j].x > result.avg.box.xmax) {
result.avg.pts[j].x = result.avg.pts[j].x - xlen;
    }
    if (result.avg.pts[j].y < result.avg.box.ymin) {
result.avg.pts[j].y = result.avg.pts[j].y + ylen;
    }
    if (result.avg.pts[j].y > result.avg.box.ymax) {
result.avg.pts[j].y = result.avg.pts[j].y - ylen;
    }
    if (result.avg.pts[j].z < result.avg.box.zmin) {
result.avg.pts[j].z = result.avg.pts[j].z + zlen;
    }
    if (result.avg.pts[j].z > result.avg.box.zmax) {
result.avg.pts[j].z = result.avg.pts[j].z - zlen;
    }
  }
  result.avg.box.xlen = xlen;
  result.avg.box.ylen = ylen;
  result.avg.box.zlen = zlen;


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

  // std::string outfile = "test.traj";
  // std::ofstream test_out{outfile};
  // for (int i =0; i < 250; i++){
  //   test_out << result.avg.pts[i].x << " " << result.avg.pts[i].y << " " << result.avg.pts[i].y << std::endl;
  // }
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


template <typename T>
void NNeighbors(std::string &filename, const size_t N, const int num_frames,
    const int num_nbs) {
  PointCloud<T> cloud;
  double skin = 4.0;
  size_t header = 5;
  int n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist,
      zdist, variance;
  double diff_sqrd = 0;
  double old_avg = 0;
  double avg = 0;
  Trajectory traj(filename, N, header);
  simFrame<T> frame;
  frame.index = 0;
  for (int i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    populatePointCloudPBC(cloud, frame, skin, true);
    typedef KDTreeSingleIndexAdaptor<
      //L2_Simple_Adaptor<T, PointCloud<T> >,
      NN_Adaptor<T, PointCloud<T> >,
      PointCloud<T>,
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
    xlen = frame.box.xlen;
    ylen = frame.box.ylen;
    zlen = frame.box.zlen;
    for (int j = 0; j < N; j++)
    {
      idx = j;
      const T query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
          cloud.pts[idx].z };
      {
      std::vector<size_t> ret_index(num_results);
      std::vector<T> out_dist(num_results);
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

template <typename T>
void variance01kd(std::string &filename, simFrame<T> &avg_frame, const size_t N,
    const int num_frames, const int num_skipframes, const int num_nbs) {
// create PointCloud object
  PointCloud<T> cloud;
  // skin defines amount of padding to add to outsides of simulation cube
  // units are Angstroms
  // Instead of the minimum image criterium which won't work for KDtree
  // we replicate enough of the region from PBC necessary to perform calculations
  double skin = 5.0;
  size_t header = 5;
  long n = 1;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist, ydist, zdist,
      variance;
  double diff_sqrd = 0.0;
  double old_avg = 0.0;
  double avg = 0.0;
  Trajectory traj(filename, N, header);
  simFrame<T> frame;
  //std::cout << "populating point cloud" << std::endl;
  // is avg_frame scaled or no?
  populatePointCloudPBC(cloud, avg_frame, skin);
  typedef KDTreeSingleIndexAdaptor<
   NN_Adaptor<T, PointCloud<T> >,
   PointCloud<T>,
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
    const T query_pt[3] = { cloud.pts[idx].x, cloud.pts[idx].y,
        cloud.pts[idx].z};
    {
    std::vector<size_t> ret_index(num_results);
    std::vector<T> out_dist(num_results);
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
    xlen = frame.box.xlen;
    ylen = frame.box.ylen;
    zlen = frame.box.zlen;
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

// template <typename T>
// NeigborList<size_t> NeighborListGenerator(PointCloud<T> cloud, size_t N,
//                     const int num_nbs,){
//
//
// }

template <typename T>
void variance01kd_r(std::string &filename, simFrame<T> &avg_frame, const size_t N,
    const int num_frames, const int num_skipframes, const int num_nbs,
    size_t &nbs_found, double &variance01, std::string &outfile, double skin,
    int dump=0, float eps=0.0001) {
  // std::cout << dump << "\n";
  // PointCloud<T> cloud;
  if (!outfile.empty()) {
    std::ofstream var_out {outfile};
    var_out << "npairs\t   variance\n";
  }

  // double skin = 6.0;
  size_t header = 5;
  size_t n = 0;
  double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist_2, ydist_2, zdist_2,
      dist, variance;
  double diff_sqrd = 0.0;
  double old_avg = 0.0;
  double avg = 0.0;
  double var_check;
  Trajectory traj(filename, N, header);
  simFrame<T> frame;

  std::vector<size_t> idxs;
  std::vector<size_t> neigh_idxs;

  NeighborListGenerator NL;
  NL = NeighborListGenerator(avg_frame, N, num_nbs, skin);
  std::vector<std::vector<size_t>> nbs = NL.GetListOG();
  // std::vector<std::vector<size_t>> nbs = nlist.nbs_og;
  // neigh_idxs = nlist.nbs;
  // size_t nbs = nlist.nnbs;
  RunningStat rs;
  traj.skipFrames(num_skipframes);
  std::string buff = "";
  int lines = 0;
  int frames_counted = 0;
  for (size_t i = 0; i < num_frames; i++) {
    traj.getNextFrame(frame);
    // std::cout << "frame number = " << i << std::endl;
    frames_counted++;
    xlen = frame.box.xlen;
    ylen = frame.box.ylen;
    zlen = frame.box.ylen;
    // std::cout << "nbs = " << nbs << std::endl;
    for (size_t j = 0; j < N; j++) {
      for (size_t k = 0; k < num_nbs; k++) {
        // grab our reference point x,y,z coordinates and coordinates of it's
        // neighbors
        xa = frame.pts[i].x;
        xb = frame.pts[nbs[i][k]].x;
        ya = frame.pts[i].y;
        yb = frame.pts[nbs[i][k]].y;
        za = frame.pts[i].z;
        zb = frame.pts[nbs[i][k]].z;
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

        // debuggin line?
        var_check = rs.Variance();

        rs.Push(dist);
        n++;
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
            if (lines < 1000) {
              std::stringstream stream1;
              stream1 << std::fixed << std::setprecision(10)
              << rs.NumDataValues();
              std::stringstream stream2;
              stream2 << std::fixed << std::setprecision(10)
              << rs.Variance();
              buff = buff + stream1.str() + " " + stream2.str() + "\n";
              lines++;
            }
            else {
              std::ofstream var_out {outfile, std::ios_base::app};
              var_out << std::fixed << std::setprecision(10) << buff;
              lines = 0;
              buff = "";
            }
          }
        }
        // variance_vector.push_back(rs.Variance());

      }
    }

  }
  if (dump > 0) {
    std::ofstream var_out {outfile, std::ios_base::app};
    var_out << std::fixed << std::setprecision(10) << buff;
  }

  nbs_found = rs.NumDataValues();
  variance01 = rs.Variance();


}




#endif
