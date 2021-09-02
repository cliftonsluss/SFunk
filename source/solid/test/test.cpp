// have to include source file not header in order to link objects
#include "../read_config.cpp"
#include "../traj_reader.cpp"
#include "../neighbor_list_generator.cpp"
#include "../test.h"
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

TEST(TrajectoryTest, GetSkipFrames) {
  std::string filename = "traj_256a_1000f_fcc.trj";
  const size_t num_atoms = 256;
  const size_t header = 5;
  Trajectory traj(filename, num_atoms, header);
  simFrame<double> frame;
  traj.getNextFrame(frame);
  double len = 7.2298*2.0;
  double min = -7.2298;
  double x = 0.9996583*len + min;
  double y = 0.9999705*len + min;
  double z = 0.0000811*len + min;
  // test that we read points correctly
  ASSERT_EQ(frame.pts[0].x, x);
  ASSERT_EQ(frame.pts[0].y, y);
  ASSERT_EQ(frame.pts[0].z, z);
  // check that we can read points correctly
  // after a skipFrames call
  traj.skipFrames(5);
  traj.getNextFrame(frame);
  x = 0.0001457*len + min;
  y = 0.9999864*len + min;
  z = 0.0001194*len + min;
  ASSERT_EQ(frame.pts[0].x, x);
  ASSERT_EQ(frame.pts[0].y, y);
  ASSERT_EQ(frame.pts[0].z, z);
}
class NeighborListTest : public ::testing::Test {
protected:
  void SetUp() override {
    num_nbs_bcc = 8;
    datafile_bcc = "traj_250a_1000f_bcc.trj";
    datafile_Fe = "Fe_dummy_2.trj";
    num_atoms_Fe = 93312;
    num_frames_Fe = 2;
    num_atoms_bcc = 250;
    num_frames = 1000;
    num_skipframes = 0;
    skin = 6.0;
    // variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
    //   results_bcc);
    // NL = NeighborListGenerator(results_bcc.avg, num_atoms_bcc,
    //                            num_nbs_bcc, skin);

  }
std::string datafile_bcc;
std::string datafile_Fe;
size_t num_atoms_bcc;
size_t num_atoms_Fe;
int num_nbs_bcc;
int num_frames;
int num_frames_Fe;
int num_skipframes;
size_t neighbors_bcc;
size_t neighbors_Fe;
resultSet<double> results_bcc;
resultSet<double> results_Fe;
double skin;
NeighborListGenerator NL_bcc;
NeighborList<size_t> nlist_bcc;
NeighborListGenerator NL_Fe;
NeighborList<size_t> nlist_Fe;
};

TEST_F(NeighborListTest, BCC_neighbors){
  variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
    results_bcc);
  NL_bcc = NeighborListGenerator(
                        results_bcc.avg, num_atoms_bcc, num_nbs_bcc,
                        skin);
  neighbors_bcc = num_atoms_bcc * num_nbs_bcc/2;
  nlist_bcc = NL_bcc.GetList();
  ASSERT_EQ(nlist_bcc.nnbs, neighbors_bcc);
}

TEST_F(NeighborListTest, Fe_neighbors){
  variance00WK<double>(datafile_Fe, num_atoms_Fe, num_frames_Fe, num_skipframes,
    results_Fe);
  NL_Fe = NeighborListGenerator(
                        results_Fe.avg, num_atoms_Fe, num_nbs_bcc,
                        skin);
  neighbors_Fe = num_atoms_Fe * num_nbs_bcc/2;
  nlist_Fe = NL_Fe.GetList();
  ASSERT_EQ(nlist_Fe.nnbs, neighbors_Fe);
}

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
    num_frames_Fe = 2;
    datafile_fcc = "traj_256a_1000f_fcc.trj";
    datafile_bcc = "traj_250a_1000f_bcc.trj";
    datafile_Fe = "Fe_dummy_2.trj";
    outfile_fcc = "fcc.var";
    outfile_bcc = "bcc.var";
    outfile_Fe = "Fe.var";
    skin = 6.0;
    skin_Fe = 2.0;
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

// test that correct number of pairs found
// test that correct variance calculated
// variance of two random independent variable is the sum of the
// variances of the individual variables
TEST_F(variance01Test, CorrectValues_bcc) {
  variance00WK<double>(datafile_bcc, num_atoms_bcc, num_frames, num_skipframes,
    results_bcc);
  neighbors_bcc = num_atoms_bcc * num_frames * num_nbs_bcc/2;
  variance01kd_r<double>(datafile_bcc, results_bcc.avg, num_atoms_bcc,
    num_frames, num_skipframes, num_nbs_bcc, nbs_found_bcc, variance01_bcc,
    outfile_bcc, skin, dump);
  ASSERT_EQ(nbs_found_bcc, neighbors_bcc);
  std::cout << "neighbors found bcc= " << nbs_found_bcc << std::endl;
  std::cout << "neighbors expected bcc= " << neighbors_bcc << std::endl;
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
  ASSERT_EQ(nbs_found_Fe, neighbors_Fe);
  std::cout << "neighbors found Fe= " << nbs_found_Fe << std::endl;
  std::cout << "neighbors expected Fe= " << neighbors_Fe << std::endl;
  ASSERT_NEAR(results_Fe.variance, variance01_Fe/2.0, 1.0e-07);
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
