#include <string>
#include <fstream>
#include <vector>
#include "sfunk_const.h"

using namespace std;

#ifndef RDF_H
#define RDF_H

namespace entropy_functional {

  class RDF {
  
    public:
      RDF(string& filename, string format="keffrdf");
      void test();
      double get_integral();
      int get_ibins() const;
      vector<double> get_r() const;
      vector<double> get_gofr() const;


    private:
      ifstream inputfile;
      vector<double> gofr;
      vector<double> r;
      vector<int> index;
      int idx;
      int ibins;
      double dr;
      int sum;
      double first_term;
      double last_term;
      double a,b,c;
      string temp;
      int integral_flag;
      double integral;
  };
}
#endif

