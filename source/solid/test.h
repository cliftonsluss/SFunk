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

#ifndef TEST_H
#define TEST_H

// using namespace std;

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
  int pbc_idx = cloud.count;
  // first we add padding to boundaries in x direction
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
  // if this point has alrady been remapped, just the map the new reference to
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
	cloud.pbc_idx_map[pbc_idx] = 1;
	pbc_idx++;
      }
    }
  }
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
    int NumDataValues() const
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
    int m_n;
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

      // the commented out section abov eand the 10 new lines below short circuit
      // original intent of this method but provide a proper Welford running avg
      // RunningStat rsx;
      // RunningStat rsy;
      // RunningStat rsz;
      //
      // rsx.Push(xb);
      // rsy.Push(yb);
      // rsz.Push(zb);
      //
      // result.avg.pts[j].x = rsx.Mean();
      // result.avg.pts[j].y = rsy.Mean();
      // result.avg.pts[j].z = rsz.Mean();
      //
      // result.variance = rsx.Variance();


      // result.variance = pow(pow(rsx.Variance(),2) + pow(rsy.Variance(),2)
          // + pow(rsz.Variance(),2), 0.5);

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




#endif
