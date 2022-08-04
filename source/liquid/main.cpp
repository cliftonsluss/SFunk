#include "funks.h"
#include "../solid/read_config.h"
#include "../solid/structures.h"
#include "rdf_reader.h"
#include <vector>
#include <iostream>

int main(int argc, char *argv[]) {

  std::string config_file = argv[1];
  std::cout << "Reading configuration file..." << std::endl;
  Read_config config(config_file);
  config.get_config();

  std::string datafile = config.j["datafile"];
  std::cout << "Reading data from " << datafile << std::endl;

  int header = config.j["header"];
  std::cout << "Skipping " << header << " header lines" << std::endl;

  double rho = config.j["density"];
  std::cout << "density = " << rho << std::endl;

  int cn = config.j["coordination_number"];
  std::cout << "coordination number = " << cn << std::endl;

  double pf = config.j["packing_fraction"];
  std::cout << "packing fraction = " << pf << std::endl;

  double c1 = config.j["c1"];
  std::cout << "c1 = " << c1 << std::endl;

  double c2 = config.j["c2"];
  std::cout << "c2 = " << c2 << std::endl;

  double q1 = config.j["q1"];
  std::cout << "q1 = " << q1 << std::endl;

  double q2 = config.j["q2"];
  std::cout << "q2 = " << q2 << std::endl;


  RDF<double> rdf;


  RDF_reader rdf_reader(datafile, header);
  rdf_reader.getData(rdf);

  SFunk sfunk(rdf, rho);

  double kirk = sfunk.kirkwood();
  double Q = sfunk.Q();
  double phit = sfunk.phi_tilde(Q, q1, q2);
  double epft = sfunk.EPFT(cn, pf, c1, c2, q1, q2);

  std::cout << "Kirkwood entropy = " << kirk << std::endl;
  std::cout << "Q = " << Q << std::endl;
  std::cout << "phi tilde = " << phit << std::endl;
  std::cout << "EPFT = " << epft << std::endl;


}
