#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include "../solid/structures.h"

#ifndef RDF_READER_H
#define RDF_READER_H



class RDF_reader {
  public:
    RDF_reader(std::string &filename,
               const int header);

    void getData(RDF<double> &rdf);

  private:
    // string filename;
    int header;
    int footer;
    std::ifstream inputfile;

};

#endif
