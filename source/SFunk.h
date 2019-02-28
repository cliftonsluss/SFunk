#include <vector>
#include "rdf.h"
#include "sfunk_const.h"
#include "sim.h"

using namespace std;

#ifndef SFUNK_H
#define SFUNK_H

namespace entropy_functional {

  class SFunk {
    public:
      SFunk(const RDF& rdf, const SIM& sim, double T);
      //double get_S(void);
      double get_S_excess_over_NkB(void);
      double get_S_excess(void);
      double get_Q(void);

    private:
      double dbw2;
      double dbw;
      double S_first;
      double S_second;
      double S_third;
      double term_first;
      double term_second;
      double term_last;
      double dr;
      double S_excess_over_NkB;
      double S_over_NkB;
      double S;
      double sum;
      double Q;
      double phi;
      double kappa_4;
      double rho;
      int N;
      double R;
      double mass;
      int ibins;
      vector<double> r;
      vector<double> gofr;
  };
}

#endif
