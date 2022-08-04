#include <cstdlib>
#include <vector>
#include <cmath>
#include "constants.h"
#include "../solid/structures.h"


#ifndef FUNKS_H
#define FUNKS_H

class SFunk {
  public:
    SFunk(RDF<double> &rdf, double rho);


    RDF<double> rdf;
    double rho;

    double kirkwood();
    double Q();
    double Kappa();
    double phi_tilde(double Q, double q1,double q2);
    double EPFT(int cn,
                double pf,
                double c1,
                double c2,
                double q1,
                double q2);

  private:
    double Prob_r(double r, double Rs);
    double Facpp(double kappa,
                 int cn,
                 int cnb,
                 double pf,
                 double pfb,
                 double rho,
                 double c1,
                 double c2);

  //   RDF<double> rdf;
  //   double rho;
};

#endif
