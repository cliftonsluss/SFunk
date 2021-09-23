#include "PBCPointCloud_generator.h"


PBCPointCloudGenerator::PBCPointCloudGenerator(simFrame<double> &frame,
                        double skin, bool scaled,
                        const double max_range){

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
  // pad outside of cloud to mimic periodic boundary conditions
  // and create map of added indexes back to original indexes
  // pbc_idx is the index of an atom at it's position in the skin
  // a map of key:value pairs is created with pbc_idx as the key and the
  // original index as the value. Any time an atom is added to the skin more
  // than once it's new pbc_idx is added to pbc_idx_map as the key and the value
  // stored at the last previous pbc_idx is stored as the value.
  int pbc_idx = cloud.count;
  // std::cout << "pbc_idx start = " << pbc_idx << std::endl;
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
      pbc_idx++;
      // std::cout << "x of index " << i << " modified by adding " << xlen
      // << " and added as new index " << pbc_idx-1 << std::endl;
    }
  }
  // std::cout << pbc_idx-diff_idx << " atoms added in x dir\n";
  diff_idx = pbc_idx;
  // we have to reset cloud.count to include the atoms we've added
  cloud.count = pbc_idx;

  // next we pad the y boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].y > y_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y - ylen,
                           cloud.pts[i].z});
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
  // std::cout << pbc_idx-diff_idx << " atoms added in y dir\n";
  diff_idx = pbc_idx;

  cloud.count = pbc_idx;
  // finally pad the z boundaries
  for (int i =0; i < cloud.count; i++) {
    if (cloud.pts[i].z > z_thresh_max) {
      cloud.pts.push_back({cloud.pts[i].x,
                           cloud.pts[i].y,
                           cloud.pts[i].z - zlen});

      if (cloud.pbc_idx_map.find(i) != cloud.pbc_idx_map.end()) {
	      cloud.pbc_idx_map[pbc_idx] = cloud.pbc_idx_map[i];
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
  // std::cout << pbc_idx-diff_idx << " atoms added in z dir\n";
  // std::cout << cloud.pbc_idx_map.size() << " atoms in padding" << std::endl;
  // std::cout << "pbc_idx = " << pbc_idx << std::endl;

};

PointCloud<double> PBCPointCloudGenerator::GetCloud(){
  return cloud;
};
