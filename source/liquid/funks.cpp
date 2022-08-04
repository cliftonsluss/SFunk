#include "funks.h"

SFunk::SFunk(RDF<double> &rdf, double rho) {

  SFunk::rdf = rdf;
  SFunk::rho = rho;
}

double SFunk::kirkwood() {
  double result;
  double tot = 0.0;
  double dr = rdf.r[1] - rdf.r[0];
  double first = 0.0;
  double last = 0.0;
  double first_term = 0.0;
  double second_term = 0.0;
  double integral;

  // first term in Kirkwood
  for(int i = 1; i < (rdf.r.size()-1); i++) {
    std::cout << rdf.r[i] << ", " << rdf.g_of_r[i] << std::endl;
    if (rdf.g_of_r[i] > 0) {
      tot = tot + rdf.g_of_r[i]*log(rdf.g_of_r[i])*pow(rdf.r[i],2);
    }
  }
  if (rdf.g_of_r.front() > 0) {
    first = rdf.g_of_r.front()*log(rdf.g_of_r.front())*pow(rdf.r.front(),2);
  }
  if (rdf.g_of_r.back() > 0) {
    first = rdf.g_of_r.back()*log(rdf.g_of_r.back())*pow(rdf.r.back(),2);
  }
  first_term = 0.5*dr*(first + 2.0*tot + last);

  // second term in kirkwood
  tot = 0.0;
  for(int i = 1; i < (rdf.r.size()-1); i++) {
    tot = tot + (rdf.g_of_r[i] -1)*pow(rdf.r[i],2);
  }
  first = (rdf.g_of_r.front() -1)*pow(rdf.r.front(),2);
  last = (rdf.g_of_r.back() -1)*pow(rdf.r.back(),2);
  second_term = 0.5*dr*(first + 2.0*tot + last);

  integral = -1.0 -0.5*(-1 + rho*4.0*PI*(first_term - second_term));
  std::cout << "pi = " << PI << std::endl;
  return integral;
}

double SFunk::Q(){
  double dr = rdf.r[1] - rdf.r[0];
  double tot = 0.0;
  double fact = rho*4.0*PI;
  double first = 0.0;
  double last = 0.0;
  double Rs = pow(0.75/PI/rho,0.3333333333);

  for(int i = 1; i < (rdf.r.size()-1); i++) {
    tot = tot + (fact*rdf.g_of_r[i]*Prob_r(rdf.r[i],Rs)*(rdf.r[i]*rdf.r[i]));
  }
  first = fact*rdf.g_of_r.front()*Prob_r(rdf.r.front(),Rs)*(rdf.r.front()*rdf.r.front());
  last = fact*rdf.g_of_r.back()*Prob_r(rdf.r.back(),Rs)*(rdf.r.back()*rdf.r.back());
  return 0.5*dr*(first + 2.0*tot + last);
}

double SFunk::Kappa() {
  double dr = rdf.r[1] - rdf.r[0];
  int even = 0;
  int odd = 0;
  double first = 0.0;
  double last = 0.0;
  double integral;

  for(int e = 1; e < (rdf.r.size()-1); e = e+2){
    even = even + pow((rdf.g_of_r[e]-1.0),2)*pow(rdf.r[e],2);
  }
  for(int o = 2; o < (rdf.r.size()-2); o = o+2){
    odd = odd + pow((rdf.g_of_r[o]-1.0),2)*pow(rdf.r[o],2);
  }

  first = pow((rdf.g_of_r.front()-1.0),2)*pow(rdf.r.front(),2);
  last = pow((rdf.g_of_r.back()-1.0),2)*pow(rdf.r.back(),2);

  integral = (dr/3.0)*rho*4.0*PI*(first + 4.0*even + 2.0*odd + last);
  return integral;
}

double SFunk::phi_tilde(double Q, double q1, double q2) {
  return Q + q1*Q*(1.0- Q) + q2*pow(Q,2)*(1.0 -Q);
}

double SFunk::EPFT(int cn,
                   double pf,
                   double c1,
                   double c2,
                   double q1,
                   double q2) {
  // we use the coordination number and packing fraction
  // for silicon as our reference CN and PF values.
  int cnb = 4;
  double pfb = 0.359619587422596;

  return kirkwood()*Facpp(Kappa(),
                          cn,
                          cnb,
                          pf,
                          pfb,
                          rho,
                          c1,
                          c2) + 0.5*phi_tilde(Q(), q1, q2);
}

double SFunk::Prob_r(double r, double Rs) {
  return exp(-(r/Rs)*(r/Rs)*(r/Rs));
}

double SFunk::Facpp(double kappa,
             int cn,
             int cnb,
             double pf,
             double pfb,
             double rho,
             double c1,
             double c2) {
  double term1;
  double term2;

  if (cn < 7 && cn > 5) {
    return 1.0;
  }
term1 = ((cnb - 6)/( cn - 6))*pow(abs((cnb - 6)/(cn - 6)),2);
term2 = pfb/pf*c1/exp(pow(c2*rho/kappa,2));
return 1.0 + term1*term2;
}
