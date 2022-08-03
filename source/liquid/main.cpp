#include "funks.h"
// #include "rdf.h"
// #include "sim.h"
#include "../solid/read_config.h"
#include "../solid/structures.h"
#include "rdf_reader.h"
#include <vector>
#include <iostream>

// using namespace std;
// using namespace entropy_functional;

int main(int argc, char *argv[]) {

  std::string config_file = argv[1];
  std::cout << "Reading configuration file..." << std::endl;
  Read_config config(config_file);
  config.get_config();

  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;

  int header = config.j["header"];
  std::cout << "Skipping " << header << " header lines" << std::endl;

  RDF<double> rdf;
  double rho = 0.43;

  RDF_reader rdf_reader(datafile, header);
  rdf_reader.getData(rdf);

  SFunk sfunk(rdf, rho);

  sfunk.kirkwood();



  // for (int i = 0; i < rdf.r.size(); i++){
  //
  //   std::cout << i << ", " << rdf.r[i] << ", " << rdf.g_of_r[i] << std::endl;
  //
  // }
  // std::cout << rdf.r[0] << ", " << rdf.g_of_r[0] << std::endl;


  // 	//double T = 1000000.0;
	// double a = 2.86997;
	// double T = 1000000.0;
	// double box = 20.0;
	// int N = 16000;
	// //double R = 28.5995;
	// double R = 28.5995;
  //       //double R = 28.60;
	// double mass = 3.36299e-21;
	// string filename = "./rdfs/iron_1000000K.rdf";
	// RDF rdf(filename);
	// rdf.test();
	// SIM sim(a, box, N, mass, R);
	// SFunk s_funk(rdf, sim, T);
  //
	// cout << "S excess  = " << s_funk.get_S_excess() << "\n";
  //       cout << "Q = " << s_funk.get_Q() << "\n";
	// cout << "rdf integral = " << rdf.get_integral() << endl;

}
