#include "/Users/clifton/googletest/googletest/include/gtest/gtest.h"

class TrajectoryTest : public ::testing::Test {
protected:
  void SetUp() override {
    datafile = "traj_54a_1f_bcc.trj";
    avgdatafile = "Fe_dummy_long.trj";
    outfile = "testout.trj";
    avgoutfile = "avgout.trj";
    outcloud = "cloud.trj";
    num_atoms = 54;
    num_avg_atoms = 93312;
    num_avg_frames = 2;
    num_skipframes = 0;
    header = 5;
    skin = 6.0;
  }
  std::string datafile, outfile, outcloud, avgdatafile, avgoutfile;
  size_t num_atoms, num_avg_atoms;
  size_t header, num_avg_frames, num_skipframes;
  simFrame<double> frame;
  resultSet<double> results;
  // simFrame<double> frame2;
  PointCloud<double> cloud;
  double skin;
};

class populatePointCloudPBCTest : public ::testing::Test {
protected:
  void SetUp() override {
    datafile_bcc = "traj_54a_1f_bcc.trj";
    num_atoms_bcc = 54;
    datafile_Fe = "Fe_dummy.trj";
    num_atoms_Fe = 93312;
    num_frames_Fe = 1;
    num_frames = 1;
    num_skipframes = 0;
    skin = 3.0;
  }

PointCloud<double> cloud;
PBCPointCloudGenerator PBCPCG;
std::string datafile_bcc;
std::string datafile_Fe;
size_t num_atoms_Fe;
size_t num_frames_Fe;
size_t num_atoms_bcc;
int num_frames;
int num_skipframes;
double skin;
simFrame<double> frame_bcc;
simFrame<double> frame_fe;

};

class NeighborListTest : public ::testing::Test {
protected:
  void SetUp() override {
    num_nbs_bcc = 8;
    num_nbs_fcc = 12;
    num_nbs_Si = 4;
    // datafile_bcc = "traj_250a_1000f_bcc.trj";
    datafile_bcc = "traj_54a_1f_bcc.trj";
    datafile_Fe = "Fe_dummy_long_2.trj";
    datafile_Si = "Si_0K.trj";
    datafile_Cu = "Cu_dummy.trj";
    num_atoms_Fe = 93312;
    num_atoms_Si = 54872;
    num_atoms_Cu = 87808;
    num_frames_Fe = 1;
    num_frames_Si = 1;
    num_frames_Cu = 1;
    // num_atoms_bcc = 250;
    num_atoms_bcc = 54;
    num_frames = 1;
    num_skipframes = 1;
    skin = 6.0;
    // variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
    //   results_bcc);
    // NL = NeighborListGenerator(results_bcc.avg, num_atoms_bcc,
    //                            num_nbs_bcc, skin);

  }
std::string datafile_bcc;
std::string datafile_Fe;
std::string datafile_Si;
std::string datafile_Cu;
size_t num_atoms_bcc;
size_t num_atoms_Fe;
size_t num_atoms_Si;
size_t num_atoms_Cu;
int num_nbs_bcc;
int num_nbs_Si;
int num_nbs_fcc;
int num_frames;
int num_frames_Fe;
int num_frames_Si;
int num_frames_Cu;
int num_skipframes;
size_t neighbors_bcc;
size_t neighbors_Fe;
size_t neighbors_Si;
size_t neighbors_Cu;
resultSet<double> results_bcc;
resultSet<double> results_Fe;
resultSet<double> results_Si;
resultSet<double> results_Cu;
double skin;
NeighborListGenerator NL_bcc;
std::vector<std::vector<size_t>> nlist_bcc;
NeighborListGenerator NL_Fe;
std::vector<std::vector<size_t>> nlist_Fe;
NeighborListGenerator NL_Si;
std::vector<std::vector<size_t>> nlist_Si;
NeighborListGenerator NL_Cu;
std::vector<std::vector<size_t>> nlist_Cu;
simFrame<double> frame_bcc;
simFrame<double> frame_fe;
simFrame<double> frame_si;
simFrame<double> frame_cu;
};


class Var01_test : public ::testing::Test {
protected:
  void SetUp() override {
    lattice_file = "avgout.trj";
    datafile = "Fe_dummy_long_2.trj";
    num_atoms = 93312;
    num_frames = 10;
    num_skipframes = 0;
    num_nbs = 8;
    nbs_found = 0;
    outfile = "runningVar.out";
    skin = 6.0;
    dump = 10;
    Trajectory traj(lattice_file, num_atoms, 5);
    traj.getNextFrame(frame);
    NeighborListGenerator NLG = NeighborListGenerator(frame,
                                num_atoms,
                                num_nbs,
                                skin);
    nbl = NLG.GetListOG();

    lattice_file_Si = "avg_Si_0.1T.trj";
    datafile_Si = "Si_0.1T_20.trj";
    num_atoms_Si = 54872;
    num_frames_Si = 20;
    num_skipframes = 0;
    num_nbs_Si = 4;
    nbs_found_Si = 0;
    outfile = "runningVar.out";
    skin_Si = 3.0;
    dump = 10;
    Trajectory traj_Si(lattice_file_Si, num_atoms_Si, 5);
    traj_Si.getNextFrame(frame_Si);

    // NeighborListGenerator NLG_Si = NeighborListGenerator(frame_Si,
    //                             num_atoms_Si,
    //                             num_nbs_Si,
    //                             skin_Si);
    // nbl_Si = NLG_Si.GetListOG();
  }
  std::string datafile;
  std::string lattice_file;
  size_t num_atoms;
  int num_frames;
  int num_skipframes;
  int num_nbs;
  size_t nbs_found;
  double variance01;
  std::vector<double> variance01_xyz;
  std::string outfile;
  double skin;
  std::vector<std::vector<size_t>> nbl;
  int dump;
  simFrame<double> frame;

  std::string datafile_Si;
  std::string lattice_file_Si;
  size_t num_atoms_Si;
  int num_frames_Si;
  int num_nbs_Si;
  size_t nbs_found_Si;
  double variance01_Si;
  std::vector<double> variance01_xyz_Si;
  double skin_Si;
  std::vector<std::vector<size_t>> nbl_Si;
  simFrame<double> frame_Si;
};

// create test fixture for variance01 tests
class variance01Test : public ::testing::Test {
protected:
  void SetUp() override {
    num_nbs_fcc = 12;
    num_nbs_bcc = 8;
    num_atoms_fcc = 256;
    num_atoms_bcc = 250;
    num_frames = 1000;
    num_skipframes = 0;
    num_atoms_Fe = 93312;
    num_frames_Fe = 10;
    datafile_fcc = "traj_256a_1000f_fcc.trj";
    datafile_bcc = "traj_250a_1000f_bcc.trj";
    datafile_Fe = "Fe_dummy_long_2.trj";
    outfile_fcc = "fcc.var";
    outfile_bcc = "bcc.var";
    outfile_Fe = "Fe.var";
    skin = 6.0;
    skin_Fe = 6.0;
    dump = 1;
    error = 0.00000001;

    // variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
    //   num_skipframes, num_nbs, nbs_found, variance01, outfile, dump);
  }
  std::string config_file;
  int num_nbs_fcc;
  int num_nbs_bcc;
  int num_atoms_fcc;
  int num_atoms_bcc;
  int num_atoms_Fe;
  int num_frames;
  int num_frames_Fe;
  int num_skipframes;
  std::string datafile_fcc;
  std::string datafile_bcc;
  std::string datafile_Fe;
  resultSet<double> results_fcc;
  resultSet<double> results_bcc;
  resultSet<double> results_Fe;
  size_t nbs_found_bcc;
  size_t nbs_found_fcc;
  size_t nbs_found_Fe;
  double variance01_fcc;
  double variance01_bcc;
  double variance01_Fe;
  size_t neighbors_fcc;
  size_t neighbors_bcc;
  size_t neighbors_Fe;
  std::vector<double> var_vec_fcc;
  std::vector<double> var_vec_bcc;
  std::vector<double> var_vec_Fe;
  std::string outfile_fcc;
  std::string outfile_bcc;
  std::string outfile_Fe;
  double skin;
  double skin_Fe;
  int dump;
  double error;
};
