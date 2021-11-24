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

// create test fixture for PBC_XYZ_Distance tests
class PBC_XYZ_Distance : public ::testing::Test {
protected:
  void SetUp() override {
    configfile_Cu = "./PBC_XYZ_Cu.json";
    Read_config config_Cu(configfile_Cu);
    config_Cu.get_config();
    num_nbs_Cu = config_Cu.j["neighbors"];
    num_atoms_Cu = config_Cu.j["atoms"];
    num_frames_Cu = config_Cu.j["frames"];
    num_skipframes_Cu = config_Cu.j["skipframes"];
    datafile_Cu = config_Cu.j["datafile"];
    skin_Cu = config_Cu.j["skin"];
    dump_Cu = config_Cu.j["dump"];
    outfile_Cu = config_Cu.j["outfile"];

    configfile_Si = "./PBC_XYZ_Si.json";
    Read_config config_Si(configfile_Si);
    config_Si.get_config();
    num_nbs_Si = config_Si.j["neighbors"];
    num_atoms_Si = config_Si.j["atoms"];
    num_frames_Si = config_Si.j["frames"];
    num_skipframes_Si = config_Si.j["skipframes"];
    datafile_Si = config_Si.j["datafile"];
    skin_Si = config_Si.j["skin"];
    dump_Si = config_Si.j["dump"];
    outfile_Si = config_Si.j["outfile"];

    // variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
    //   num_skipframes, num_nbs, nbs_found, variance01, outfile, dump);
  }
std::string config_file;
int num_nbs_Cu;
int num_atoms_Cu;
int num_frames_Cu;
int num_skipframes_Cu;
std::string datafile_Cu;
std::string configfile_Cu;
resultSet<double> results_Cu;
size_t nbs_found_Cu;
double variance01_Cu;
size_t neighbors_Cu;
std::string outfile_Cu;
double skin_Cu;
int dump_Cu;
std::vector<std::vector<size_t>> nbl_Cu;

int num_nbs_Si;
int num_atoms_Si;
int num_frames_Si;
int num_skipframes_Si;
std::string datafile_Si;
std::string configfile_Si;
resultSet<double> results_Si;
size_t nbs_found_Si;
double variance01_Si;
size_t neighbors_Si;
std::string outfile_Si;
double skin_Si;
int dump_Si;
std::vector<std::vector<size_t>> nbl_Si;

std::vector<double> variance_xyz;
NeighborListGenerator NL;
double xlen, ylen, zlen, xa, ya, za, xb, yb, zb, xdist_2, ydist_2, zdist_2,
    dist, variance;
std::vector<double> dist_xyz;
RunningStat rs;
RunningStat rsx;
RunningStat rsy;
RunningStat rsz;
};


// create test fixture for variance01 tests
class variance01Test : public ::testing::Test {
protected:
  void SetUp() override {
    configfile_Cu = "./Cu.json";
    Read_config config_Cu(configfile_Cu);
    config_Cu.get_config();
    num_nbs_Cu = config_Cu.j["neighbors"];
    num_atoms_Cu = config_Cu.j["atoms"];
    num_frames_Cu = config_Cu.j["frames"];
    datafile_Cu = config_Cu.j["datafile"];
    skin_Cu = config_Cu.j["skin"];
    dump_Cu = config_Cu.j["dump"];
    outfile_Cu = config_Cu.j["outfile"];


    num_nbs_fcc = 12;
    num_nbs_bcc = 8;
    num_nbs_Si = 4;
    num_nbs_Fe = 4;
    num_atoms_fcc = 256;

    num_atoms_bcc = 250;
    num_frames = 1000;
    num_skipframes = 0;
    num_atoms_Fe = 93312;
    num_atoms_Si = 54872;
    num_frames_Fe = 10;
    num_frames_Si = 20;

    datafile_fcc = "traj_256a_1000f_fcc.trj";
    datafile_bcc = "traj_250a_1000f_bcc.trj";
    datafile_Fe = "Fe_dummy_long_2.trj";
    datafile_Si = "Si_0.8435_20.trj";

    outfile_fcc = "fcc.var";
    outfile_bcc = "bcc.var";
    outfile_Fe = "Fe.var";
    outfile_Si = "Si.var";
    // outfile_Cu = "Cu.var";
    skin = 6.0;
    skin_Fe = 6.0;
    skin_Si = 6.0;

    dump = 1;
    error = 0.00000001;

    // variance01kd_r<double>(datafile, results.avg, num_atoms, num_frames,
    //   num_skipframes, num_nbs, nbs_found, variance01, outfile, dump);
  }
  std::string config_file;
  int num_nbs_fcc;
  int num_nbs_bcc;
  int num_nbs_Si;
  int num_nbs_Cu;
  int num_nbs_Fe;
  int num_atoms_fcc;
  int num_atoms_bcc;
  int num_atoms_Fe;
  int num_atoms_Si;
  int num_atoms_Cu;
  int num_frames;
  int num_frames_Fe;
  int num_frames_Si;
  int num_frames_Cu;
  int num_skipframes;
  std::string datafile_fcc;
  std::string datafile_bcc;
  std::string datafile_Fe;
  std::string datafile_Si;
  std::string datafile_Cu;
  std::string configfile_Cu;
  resultSet<double> results_fcc;
  resultSet<double> results_bcc;
  resultSet<double> results_Fe;
  resultSet<double> results_Si;
  resultSet<double> results_Cu;
  size_t nbs_found_bcc;
  size_t nbs_found_fcc;
  size_t nbs_found_Fe;
  size_t nbs_found_Si;
  size_t nbs_found_Cu;
  double variance01_fcc;
  double variance01_bcc;
  double variance01_Fe;
  double variance01_Si;
  double variance01_Cu;
  size_t neighbors_fcc;
  size_t neighbors_bcc;
  size_t neighbors_Fe;
  size_t neighbors_Si;
  size_t neighbors_Cu;
  std::vector<double> variance_xyz;
  // std::vector<double> variance01_xyz_bcc;
  // std::vector<double> variance01_xyz_Fe;
  // std::vector<double> variance01_xyz_Si;
  // std::vector<double> variance00_xyz_Si;
  std::string outfile_fcc;
  std::string outfile_bcc;
  std::string outfile_Fe;
  std::string outfile_Si;
  std::string outfile_Cu;
  double skin;
  double skin_Fe;
  double skin_Si;
  double skin_Cu;
  int dump;
  int dump_Cu;
  double error;
  std::vector<std::vector<size_t>> nbl_Si;
  std::vector<std::vector<size_t>> nbl_Cu;
};
