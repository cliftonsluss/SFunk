#include "rdf_reader.h"

RDF_reader::RDF_reader(std::string &filename,
                       const int header) {

  RDF_reader::header = header;
  RDF_reader::inputfile.open(filename, std::ifstream::in);
}

void RDF_reader::getData(RDF<double> &rdf) {
  std::string temp;
  double r, g;
  for (size_t i = 0; i < header; i++) {
    getline(inputfile, temp);
  }

  while (inputfile >> r >> g){
    rdf.r.push_back(r);
    rdf.g_of_r.push_back(g);
  }
}
