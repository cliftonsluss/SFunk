#include "rdf.h"
#include <string>
#include <fstream>
#include <cmath>
#include <vector>


using namespace std;

namespace entropy_functional {

  RDF::RDF(string& filename, string format) {
    ifstream inputfile(filename);
  
    for (int i = 0; i < 4; i++) {
      getline(inputfile, temp);
    }

    ibins = 0; 


    if (format.compare("keffrdf") == 0) {
      while (inputfile >> a >> b) {
	r.push_back(a);
	gofr.push_back(b);
	ibins++;
      }
    }
      else {
	while (inputfile >> idx >> a >> b >> c) {
	  r.push_back(a);
	  gofr.push_back(b);
	  ibins++;
	}
      }
    }

  int RDF::get_ibins() const{
    return ibins;
  }

  vector<double>  RDF::get_r() const {
    return r;
  }

  vector<double> RDF::get_gofr() const {
    return gofr;
  }

  void RDF::test() {
    RDF::integral_flag = 1;
    dr = r[1] - r[0];
    sum = 0;
    for (int i = 1; i<ibins; i++) {
      sum = sum + gofr[i]*pow(r[i],2);
    }
    first_term = gofr[0]*pow(r[0],2);
    last_term = gofr[ibins]*pow(r[ibins],2);
    RDF::integral = 2.0*pi*dr*(first_term + 2.0*sum + last_term);
    return;
  }

  double RDF::get_integral() {
    if (RDF::integral_flag == 1) {
      return RDF::integral;
    } else {
      return -1.0;
    }
  }
}






