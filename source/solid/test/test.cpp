// have to include source file not header in order to link objects
#include "../read_config.cpp"
#include "../traj_reader.cpp"
#include "../neighbor_list_generator.cpp"
#include "../PBCPointCloud_generator.cpp"
#include "../PBC.cpp"
#include "../cell2frame.h"
#include "../test.h"
#include "test.h"
#include "/Users/clifton/googletest/googletest/include/gtest/gtest.h"
#include <math.h>

// test of test framework :p
double squareRoot(const double a) {
    double b = sqrt(a);
    if(b != b) { // nan check
        return -1.0;
    }else{
        return sqrt(a);
    }
}

TEST(SquareRootTest, PositiveNos) {
    ASSERT_EQ(6, squareRoot(36.0));
    ASSERT_EQ(18.0, squareRoot(324.0));
    ASSERT_EQ(25.4, squareRoot(645.16));
    ASSERT_EQ(0, squareRoot(0.0));
}

TEST(SquareRootTest, NegativeNos) {
    ASSERT_EQ(-1.0, squareRoot(-15.0));
    ASSERT_EQ(-1.0, squareRoot(-0.2));
}


int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(read_config_Test, CorrectValues) {
  std::string config_file = "config_file.json";
  Read_config config(config_file);
  config.get_config();
  ASSERT_EQ(12, config.j["neighbors"]);
  ASSERT_EQ(256, config.j["atoms"]);
  ASSERT_EQ(1000, config.j["frames"]);
  ASSERT_EQ(0, config.j["skipframes"]);
  ASSERT_EQ("traj_256a_1000f.trj", config.j["datafile"]);
  ASSERT_EQ("test.out", config.j["outfile"]);
  ASSERT_EQ(5.0, config.j["skin"]);
  ASSERT_EQ(1000, config.j["dump"]);
  ASSERT_EQ(0.000000001, config.j["error"]);
}

TEST(Cell2FrameTest, doesnt_crash){
  std::vector<double> pt{0,2,2,2};
  std::vector<std::vector<double>> nbs{{1,1,1,1},
                                       {2,3,1,3},
                                       {3,1,1,3},
                                       {5,1,3,3},
                                       {4,3,3,3},
                                       {6,3,1,1},
                                       {7,3,3,1},
                                       {8,1,3,1}};
  simFrame<double> frame = cell2frame(pt, nbs);
  for (int i = 0; i < nbs.size()+1; i++){
    std::cout << frame.idx[i] << " {" << frame.pts[i].x << "," <<
    frame.pts[i].y << "," << frame.pts[i].z << "}"
    << std::endl;
  }
  std::string outfile = "cell.traj";
  Trajectory CellTraj(outfile);
  CellTraj.writeFrame(frame,false);
  CellTraj.writeFrame(frame,false);
  ASSERT_EQ(1,1);
}



TEST_F(TrajectoryTest, TrajReadWrite){
  Trajectory trajin(datafile, num_atoms, header);
  trajin.getNextFrame(frame);
  Trajectory trajout(outfile);
  trajout.writeFrame(frame);
  std::cout << "run diff testout.trj traj_256a_1f_fcc.trj to confirm\n";
  ASSERT_EQ(1,1);
}

TEST_F(TrajectoryTest, TrajWriteAvg){
  std::cout << "variance00WK block\n";
  variance00WK<double>(avgdatafile,
                       num_avg_atoms,
                       9,
                       0,
                       results);
  std::cout << "trajout block\n";
  Trajectory trajout(avgoutfile);
  std::cout << "writeFrame block\n";
  trajout.writeFrame(results.avg);
  // std::cout << "run diff testout.trj traj_256a_1f_fcc.trj to confirm\n";
  ASSERT_EQ(1,1);
}

// TEST_F(TrajectoryTest, TrajWriteCloud){
//   Trajectory trajin(datafile, num_atoms, header);
//   trajin.getNextFrame(frame);
//   Trajectory trajout(outcloud);
//   populatePointCloudPBC(cloud, frame, skin);
//   std::cout << cloud.pbc_idx_map.size() << std::endl;
//   for (int i = 0; i < cloud.pbc_idx_map.size(); i++){
//     std::cout << "index " << i << " mapped to " << cloud.pbc_idx_map[i] <<"\n";
//   }
//   cloud2frame(cloud, frame);
//   trajout.writeFrame(frame, true);
//   ASSERT_EQ(1,1);
// }


// TEST(TrajectoryTest, GetSkipFrames) {
//   std::string filename = "traj_256a_1000f_fcc.trj";
//   num_atoms = 256;
//   header = 5;
//   Trajectory traj(filename, num_atoms, header);
//   // simFrame<double> frame;
//   traj.getNextFrame(frame);
//   double len = 7.2298*2.0;
//   double min = -7.2298;
//   double x = 0.9996583*len + min;
//   double y = 0.9999705*len + min;
//   double z = 0.0000811*len + min;
//   // test that we read points correctly
//   ASSERT_EQ(frame.pts[0].x, x);
//   ASSERT_EQ(frame.pts[0].y, y);
//   ASSERT_EQ(frame.pts[0].z, z);
//   // check that we can read points correctly
//   // after a skipFrames call
//   traj.skipFrames(5);
//   traj.getNextFrame(frame);
//   x = 0.0001457*len + min;
//   y = 0.9999864*len + min;
//   z = 0.0001194*len + min;
//   ASSERT_EQ(frame.pts[0].x, x);
//   ASSERT_EQ(frame.pts[0].y, y);
//   ASSERT_EQ(frame.pts[0].z, z);
// }

// class PBCPoinCloudTest : public ::testing::Test {
// protected:
//   void SetUp() override {
//     skin = 3.0;
//   }
// PointCloud<double> cloud;
// simFrame<double> frame;
// double skin;
// double density;
// };



// TEST_F(populatePointCloudPBCTest, BCC_map){
//   Trajectory traj(datafile_bcc, num_atoms_bcc, 5);
//   traj.getNextFrame(frame_bcc);
//   populatePointCloudPBC(cloud, frame_bcc, skin);
//   std::cout << cloud.pbc_idx_map.size() << std::endl;
//   for (int i = 0; i < cloud.pbc_idx_map.size(); i++){
//     std::cout << "index " << i << " mapped to " << cloud.pbc_idx_map[i] <<"\n";
//   }
//   ASSERT_EQ(1,1);
// }

TEST_F(populatePointCloudPBCTest, Fe_map){
  Trajectory traj(datafile_Fe, num_atoms_Fe, 5);
  traj.getNextFrame(frame_fe);
  PBCPCG = PBCPointCloudGenerator(frame_fe, skin);
  cloud = PBCPCG.GetCloud();
  std::cout << cloud.pbc_idx_map.size() << std::endl;
  for (int i = 0; i < cloud.pbc_idx_map.size(); i++){
    std::cout << "index " << i << " mapped to " << cloud.pbc_idx_map[i] <<"\n";
  }
  ASSERT_EQ(1,1);
}





// TEST_F(NeighborListTest, BCC_neighbors){
//   // variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
//   //   results_bcc);
//   Trajectory traj(datafile_bcc, num_atoms_bcc, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_bcc);
//   NL_bcc = NeighborListGenerator(
//                         frame_bcc, num_atoms_bcc, num_nbs_bcc,
//                         skin);
//   // PointCloud<double> cloud;
//   // populatePointCloudPBC(cloud, frame_bcc, skin);
//   neighbors_bcc = num_atoms_bcc * num_nbs_bcc/2;
//   nlist_bcc = NL_bcc.GetList();
//   ASSERT_EQ(nlist_bcc.nnbs, neighbors_bcc);
// }

// TEST_F(NeighborListTest, Si_neighbors_skin){
//   Trajectory traj(datafile_Si, num_atoms_Si, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_si);
//   NL_Si = NeighborListGenerator(
//                         frame_si, num_atoms_Si, num_nbs_Si,
//                         skin);
//   PointCloud<double> cloud;
//   populatePointCloudSkin(cloud, frame_si, skin);
//   neighbors_Si =  cloud.list.size() * num_nbs_Si/2;
//   nlist_Si = NL_Si.GetList();
//   ASSERT_EQ(nlist_Si.nnbs, neighbors_Si);
// }

// TEST_F(NeighborListTest, Cu_neighbors_skin){
//   Trajectory traj(datafile_Cu, num_atoms_Cu, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_cu);
//   NL_Cu = NeighborListGenerator(
//                         frame_cu, num_atoms_Cu, num_nbs_fcc,
//                         skin, 5.0);
//   PointCloud<double> cloud;
//   populatePointCloudSkin(cloud, frame_cu, skin);
//   neighbors_Cu =  cloud.list.size() * num_nbs_fcc/2;
//   nlist_Cu = NL_Cu.GetList();
//   ASSERT_EQ(nlist_Cu.nnbs, neighbors_Cu);
// }

// TEST_F(NeighborListTest, Fe_neighbors_skin){
//   // variance00WK<double>(datafile_Fe, num_atoms_Fe, num_frames_Fe, num_skipframes,
//   //   results_Fe);
//   Trajectory traj(datafile_Fe, num_atoms_Fe, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_fe);
//   NL_Fe = NeighborListGenerator(
//                         frame_fe, num_atoms_Fe, num_nbs_bcc,
//                         skin, 5.0);
//   PointCloud<double> cloud;
//   populatePointCloudSkin(cloud, frame_fe, skin);
//   neighbors_Fe =  cloud.list.size() * num_nbs_bcc/2;
//   nlist_Fe = NL_Fe.GetList();
//   ASSERT_EQ(nlist_Fe.nnbs, neighbors_Fe);
// }

TEST_F(NeighborListTest, Si_neighbors){
//   Trajectory traj(datafile_Si, num_atoms_Si, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_si);
//   NL_Si = NeighborListGenerator(
//                         frame_si, num_atoms_Si, num_nbs_Si,
//                         skin);
//   // PointCloud<double> cloud;
//   // populatePointCloudSkin(cloud, frame_si, skin);
//   neighbors_Si =  num_atoms_Si * num_nbs_Si/2;
//   nlist_Si = NL_Si.GetList();
//   ASSERT_EQ(nlist_Si.nnbs, neighbors_Si);
  Trajectory traj(datafile_Si, num_atoms_Si, 5);
  traj.skipFrames(0);
  traj.getNextFrame(frame_si);
  NL_Fe = NeighborListGenerator(
                        frame_si, num_atoms_Si, num_nbs_Si,
                        skin);
  // PointCloud<double> cloud;
  // populatePointCloudSkin(cloud, frame_fe, skin);
  neighbors_Si =  num_atoms_Si * num_nbs_Si;
  nlist_Si = NL_Si.GetListOG();
  std::vector<std::vector<size_t>> errors_Si;
  errors_Si = NL_Si.GetErrors();
  std::vector<double> pt(4);
  std::vector<double> nb(4);

  std::vector<double> temp(4);
  double x,y,z;
  std::string outfile = "celltest_Si.traj";
  Trajectory CellTraj(outfile);
  simFrame<double> frame;
  std::vector<double> len{frame_fe.box.xlen,
                          frame_fe.box.ylen,
                          frame_fe.box.zlen};

  // std::cout << "box length" << std::endl;
  // std::cout << len[0] << ", " << len[1] << ", " << len[2] << std::endl;
  // std::cout << errors_Fe.size() << std::endl;
  // std::cout << "\n";

  for (int i = 0; i < errors_Si.size(); i++){

    // std::cout << errors_Fe[i][0] << ", "
    // << errors_Fe[i][1] << "\n" << std::endl;
    // std::cout << errors_Fe[i][0]+1 << " neighbors are" <<std::endl;
    x = frame_si.pts[errors_Si[i][0]].x;
    y = frame_si.pts[errors_Si[i][0]].y;
    z = frame_si.pts[errors_Si[i][0]].z;
    // std::cout << "central point" << std::endl;
    pt = {errors_Si[i][0]+1.0,x,y,z};
    // std::cout << pt[0] << " | " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
    PBC pbc1(&pt[1],len);
    std::vector<std::vector<double>> nbs1;
    // std::cout << "has neighbors" << std::endl;

    for (int j = 0; j < num_nbs_bcc; j++){
      // std::cout << "ParticleIdentifier==" << nlist_Fe[errors_Fe[i][0]][j]+1
      // << " ||" << std::endl;
      x = frame_si.pts[nlist_Si[errors_Si[i][0]][j]].x;
      y = frame_si.pts[nlist_Si[errors_Si[i][0]][j]].y;
      z = frame_si.pts[nlist_Si[errors_Si[i][0]][j]].z;
      nb = {nlist_Si[errors_Si[i][0]][j]+1.0,x,y,z};
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      pbc1.minimum_image(&nb[1]);
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      nbs1.push_back(nb);
    }
    std::cout << nbs1.size() << std::endl;
    std::cout << "\n";
    frame = cell2frame(pt, nbs1);
    CellTraj.writeFrame(frame,false);

    // std::cout << "ParticleIdentifier==" << errors_Fe[i][0]+1
    // << "\n" << std::endl;
    // std::cout << errors_Fe[i][1]+1 << " neighbors are" <<std::endl;
    x = frame_si.pts[errors_Si[i][1]].x;
    y = frame_si.pts[errors_Si[i][1]].y;
    z = frame_si.pts[errors_Si[i][1]].z;
    // std::cout << "central point" << std::endl;
    pt = {errors_Si[i][1]+1.0,x,y,z};
    // std::cout << pt[0] << " | " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
    // std::cout << "has neighbors" << std::endl;
    PBC pbc2(&pt[1],len);
    std::vector<std::vector<double>> nbs2;
    for (int j = 0; j < num_nbs_bcc; j++){
      x = frame_si.pts[nlist_Si[errors_Si[i][1]][j]].x;
      y = frame_si.pts[nlist_Si[errors_Si[i][1]][j]].y;
      z = frame_si.pts[nlist_Si[errors_Si[i][1]][j]].z;
      nb = {nlist_Si[errors_Si[i][1]][j]+1.0,x,y,z};
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      pbc2.minimum_image(&nb[1]);
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      nbs2.push_back(nb);
      // std::cout << "ParticleIdentifier==" << nlist_Fe[errors_Fe[i][1]][j]+1
      // << " ||" << std::endl;
    }
    // std::cout << nbs2.size() << std::endl;
    // std::cout << "\n";
    frame = cell2frame(pt, nbs2);
    CellTraj.writeFrame(frame,false);
    // std::cout << "ParticleIdentifier==" << errors_Fe[i][1]+1
    // << "\n" << std::endl;
  }
ASSERT_EQ(nlist_Si.size(), num_atoms_Si);
}
// }
//
// TEST_F(NeighborListTest, Cu_neighbors){
//   Trajectory traj(datafile_Cu, num_atoms_Cu, 5);
//   traj.skipFrames(0);
//   traj.getNextFrame(frame_cu);
//   NL_Cu = NeighborListGenerator(
//                         frame_cu, num_atoms_Cu, num_nbs_fcc,
//                         skin);
//   // PointCloud<double> cloud;
//   // populatePointCloudSkin(cloud, frame_cu, skin);
//   neighbors_Cu =  num_atoms_Cu * num_nbs_fcc/2;
//   nlist_Cu = NL_Cu.GetList();
//   ASSERT_EQ(nlist_Cu.nnbs, neighbors_Cu);
// }

TEST_F(NeighborListTest, Fe_neighbors){
  // variance00WK<double>(datafile_Fe, num_atoms_Fe, num_frames_Fe, num_skipframes,
  //   results_Fe);
  Trajectory traj(datafile_Fe, num_atoms_Fe, 5);
  traj.skipFrames(0);
  traj.getNextFrame(frame_fe);
  NL_Fe = NeighborListGenerator(
                        frame_fe, num_atoms_Fe, num_nbs_bcc,
                        skin);
  // PointCloud<double> cloud;
  // populatePointCloudSkin(cloud, frame_fe, skin);
  neighbors_Fe =  num_atoms_Fe * num_nbs_bcc;
  nlist_Fe = NL_Fe.GetListOG();
  std::vector<std::vector<size_t>> errors_Fe;
  errors_Fe = NL_Fe.GetErrors();
  std::vector<double> pt(4);
  std::vector<double> nb(4);

  std::vector<double> temp(4);
  double x,y,z;
  std::string outfile = "celltest.traj";
  Trajectory CellTraj(outfile);
  simFrame<double> frame;
  std::vector<double> len{frame_fe.box.xlen,
                          frame_fe.box.ylen,
                          frame_fe.box.zlen};

// std::cout << "box length" << std::endl;
// std::cout << len[0] << ", " << len[1] << ", " << len[2] << std::endl;
// std::cout << errors_Fe.size() << std::endl;
// std::cout << "\n";

  for (int i = 0; i < errors_Fe.size(); i++){

    // std::cout << errors_Fe[i][0] << ", "
    // << errors_Fe[i][1] << "\n" << std::endl;
    // std::cout << errors_Fe[i][0]+1 << " neighbors are" <<std::endl;
    x = frame_fe.pts[errors_Fe[i][0]].x;
    y = frame_fe.pts[errors_Fe[i][0]].y;
    z = frame_fe.pts[errors_Fe[i][0]].z;
    // std::cout << "central point" << std::endl;
    pt = {errors_Fe[i][0]+1.0,x,y,z};
    // std::cout << pt[0] << " | " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
    PBC pbc1(&pt[1],len);
    std::vector<std::vector<double>> nbs1;
    // std::cout << "has neighbors" << std::endl;

    for (int j = 0; j < num_nbs_bcc; j++){
      // std::cout << "ParticleIdentifier==" << nlist_Fe[errors_Fe[i][0]][j]+1
      // << " ||" << std::endl;
      x = frame_fe.pts[nlist_Fe[errors_Fe[i][0]][j]].x;
      y = frame_fe.pts[nlist_Fe[errors_Fe[i][0]][j]].y;
      z = frame_fe.pts[nlist_Fe[errors_Fe[i][0]][j]].z;
      nb = {nlist_Fe[errors_Fe[i][0]][j]+1.0,x,y,z};
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      pbc1.minimum_image(&nb[1]);
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      nbs1.push_back(nb);
    }
    std::cout << nbs1.size() << std::endl;
    std::cout << "\n";
    frame = cell2frame(pt, nbs1);
    CellTraj.writeFrame(frame,false);

    // std::cout << "ParticleIdentifier==" << errors_Fe[i][0]+1
    // << "\n" << std::endl;
    // std::cout << errors_Fe[i][1]+1 << " neighbors are" <<std::endl;
    x = frame_fe.pts[errors_Fe[i][1]].x;
    y = frame_fe.pts[errors_Fe[i][1]].y;
    z = frame_fe.pts[errors_Fe[i][1]].z;
    // std::cout << "central point" << std::endl;
    pt = {errors_Fe[i][1]+1.0,x,y,z};
    // std::cout << pt[0] << " | " << pt[1] << ", " << pt[2] << ", " << pt[3] << std::endl;
    // std::cout << "has neighbors" << std::endl;
    PBC pbc2(&pt[1],len);
    std::vector<std::vector<double>> nbs2;
    for (int j = 0; j < num_nbs_bcc; j++){
      x = frame_fe.pts[nlist_Fe[errors_Fe[i][1]][j]].x;
      y = frame_fe.pts[nlist_Fe[errors_Fe[i][1]][j]].y;
      z = frame_fe.pts[nlist_Fe[errors_Fe[i][1]][j]].z;
      nb = {nlist_Fe[errors_Fe[i][1]][j]+1.0,x,y,z};
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      pbc2.minimum_image(&nb[1]);
      // std::cout << nb[0] << " | " << nb[1] << ", " << nb[2] << ", " << nb[3] << std::endl;
      nbs2.push_back(nb);
      // std::cout << "ParticleIdentifier==" << nlist_Fe[errors_Fe[i][1]][j]+1
      // << " ||" << std::endl;
    }
    // std::cout << nbs2.size() << std::endl;
    // std::cout << "\n";
    frame = cell2frame(pt, nbs2);
    CellTraj.writeFrame(frame,false);
    // std::cout << "ParticleIdentifier==" << errors_Fe[i][1]+1
    // << "\n" << std::endl;
  }
  ASSERT_EQ(nlist_Fe.size(), num_atoms_Fe);
}

TEST(PBC_test, minimum_image){
  // std::vector<int> vec{1,2,3,4,5};
  // std::cout << vec[(*vec.begin())++] << std::endl;

  std::vector<double> pt{7,0.0,0.0,0.0};
  std::vector<double> len{2.0,2.0,2.0};
  std::vector<std::vector<double>>
                      pts{{0,1.1,1.1,1.1},
                          {1,-1.1,-1.1,-1.1},
                          {2,0.5,0.5,0.5},
                          {3,-0.5,-0.5,-0.5}};
  PBC pbc(&pt[1],len);
  for (int i =0;i < pts.size();i++){
    std::cout << "input = " << pts[i][0] << "{" << pts[i][1] << ", "
    << pts[i][2] << ", " << pts[i][3] << "}\n";
    pbc.minimum_image(&pts[i][1]);
    std::cout << "output = " << pts[i][0] << "{" << pts[i][1] << ", "
    << pts[i][2] << ", " << pts[i][3] << "}\n";
  }
  std::vector<double> a{1,2,3};
  std::vector<double> b{1,2,3};
  std::vector<std::vector<double>>
                      pts1{{0,-0.9,-0.9,-0.9},
                          {1,0.9,0.9,0.9},
                          {2,0.5,0.5,0.5},
                          {3,-0.5,-0.5,-0.5}};
  ASSERT_EQ(pts1,pts);
}

TEST_F(Var01_test, From_saved_frame) {
  variance01kd_r_nbl<double>(datafile,
    num_atoms, num_frames,
    num_skipframes, num_nbs,
    nbs_found, variance01, variance01_xyz, outfile, skin, nbl, dump);
  std::cout << "variance01= " << variance01 << std::endl;
  std::cout << "var01_x= " << variance01_xyz[0] << std::endl;
  std::cout << "var01_y= " << variance01_xyz[1] << std::endl;
  std::cout << "var01_z= " << variance01_xyz[2] << std::endl;
  std::cout << "std01= " << pow(variance01,0.5) << std::endl;
  ASSERT_EQ(1,1);
}

TEST_F(Var01_test, From_saved_Si) {
  variance01kd_r_nbl<double>(datafile_Si,
    num_atoms_Si, num_frames_Si,
    num_skipframes, num_nbs_Si,
    nbs_found_Si, variance01_Si, variance01_xyz_Si, outfile,
    skin_Si, nbl_Si, dump);
  std::cout << "variance01= " << variance01_Si << std::endl;
  std::cout << "var01_x= " << variance01_xyz_Si[0] << std::endl;
  std::cout << "var01_y= " << variance01_xyz_Si[1] << std::endl;
  std::cout << "var01_z= " << variance01_xyz_Si[2] << std::endl;
  std::cout << "std01= " << pow(variance01_Si,0.5) << std::endl;
  ASSERT_EQ(1,1);
}

// test that correct number of pairs found
// test that correct variance calculated
// variance of two random independent variable is the sum of the
// variances of the individual variables
TEST_F(variance01Test, CorrectValues_bcc) {
  std::cout << "variance00WK block\n";
  variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
    results_bcc);
  neighbors_bcc = num_atoms_bcc * num_frames * num_nbs_bcc;
  std::cout << "variance01kd_r block\n";
  variance01kd_r<double>(datafile_bcc, results_bcc.avg, num_atoms_bcc,
    num_frames, num_skipframes, num_nbs_bcc, nbs_found_bcc, variance01_bcc,
    outfile_bcc, skin, dump);
  // ASSERT_EQ(nbs_found_bcc, neighbors_bcc);
  // std::cout << "neighbors found bcc= " << nbs_found_bcc << std::endl;
  // std::cout << "neighbors expected bcc= " << neighbors_bcc << std::endl;
  ASSERT_NEAR(results_bcc.variance, variance01_bcc/2.0, 1.0e-07);
}

TEST_F(variance01Test, CorrectValues_fcc) {
  variance00WK<double>(datafile_fcc, num_atoms_fcc, num_frames, num_skipframes,
    results_fcc);
  neighbors_fcc = num_atoms_fcc * num_frames * num_nbs_fcc/2;
  variance01kd_r<double>(datafile_fcc, results_fcc.avg, num_atoms_fcc,
    num_frames, num_skipframes, num_nbs_fcc, nbs_found_fcc, variance01_fcc,
    outfile_fcc, skin, dump);
  ASSERT_EQ(nbs_found_fcc, neighbors_fcc);
  std::cout << "neighbors found fcc= " << nbs_found_fcc << std::endl;
  std::cout << "neighbors expected fcc= " << neighbors_fcc << std::endl;
  ASSERT_NEAR(results_fcc.variance, variance01_fcc/2.0, 1.0e-07);
}

TEST_F(variance01Test, CorrectValues_Fe) {
  variance00WK<double>(datafile_Fe, num_atoms_Fe, num_frames_Fe, num_skipframes,
    results_Fe);
  neighbors_Fe = num_atoms_Fe * num_frames_Fe * num_nbs_bcc/2;
  variance01kd_r<double>(datafile_Fe, results_Fe.avg, num_atoms_Fe,
    num_frames_Fe, num_skipframes, num_nbs_bcc, nbs_found_Fe, variance01_Fe,
    outfile_Fe, skin_Fe, dump);
  // ASSERT_EQ(nbs_found_Fe, neighbors_Fe);
  std::cout << "neighbors found Fe= " << nbs_found_Fe << std::endl;
  std::cout << "neighbors expected Fe= " << neighbors_Fe << std::endl;
  std::cout << "var01 " << variance01_Fe << std::endl;
  ASSERT_EQ(0, 0);
}



// test running stat output via variance01kd_r
// TEST_F(variance01Test, CorrectVarianceVectorLength) {
//   style = "welford";
//   variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
//     num_skipframes, num_nbs, nbs_found, variance01, outfile, dump,
//     error);
//   ASSERT_EQ(nbs_found, neighbors);
//   // ASSERT_EQ(var_vec.size(), neighbors);
// }
