#include "SFunk.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;

namespace entropy_functional {

  // kappa from equations 24
  // kappa = 4pi*rho*integral r^2[g(r) - 1]^2 dr  
  double kappa(double rho, vector<double> r, vector<double> gofr, int ibins, double dr) {
    double sum = 0.0;
    
    for (int i = 1; i<ibins-1; i++){
      //cout << r[i] << " " << gofr[i] << endl;
      sum = sum + pow((gofr[i]-1.0),2)*pow(r[i],2);
    }

    double term_first = pow((gofr[0]-1.0),2)*pow(r[0],2);
    double term_last  = pow((gofr[ibins-1]-1.0),2)*pow(r[ibins-1],2);
    double integral = 0.5*dr*rho*4.0*pi*(term_first + 2.0*sum + term_last);
    cout << "kappa " << integral << endl;
    return integral; 
  } 

  
  //Q from equation 29
  // N(1 -rho*integral dr(1-g_s)/(N - rho*(4/3)*pi*R^3)
  //
  double Qofg(double rho, int N, vector<double> r, vector<double> gofr, int ibins, double dr, double R) {
    double sum = 0.0;
    cout << "N inside Qofg " << N << endl;
    cout << "R inside Qofg " << R << endl;
    double denominator = 1.0 - (pow(R,3.0)*rho*(4.0/3.0)*pi/N);
    cout << "pi = " << pi << endl;
    cout << "denominator = " << denominator << endl;
    double fact = rho*4.0*pi;
   
    for (int i = 1; i<ibins -1; i++) {
      sum = sum + (1.0 - gofr[i])*pow(r[i],2);
    }

    double term_first = (1.0 - gofr[0])*pow(r[0],2);
    double term_last = (1.0 - gofr[ibins-1])*pow(r[ibins-1],2);
    double integral = 0.5*dr*(term_first + 2.0*sum + term_last);

    cout << "fact " << fact << endl;
    cout << "integra " << integral << endl;
    cout << "rho " << rho << endl;
    //return 0.5*dr*c_fac*(term_first + 2.0*sum + term_last);

    return (1.0 - fact*integral)/denominator;
  }

  // phi~ from equation 30
  double phi_tilde(double Q, double q1, double q2) {
    double T = 4000;
    double mu = 0.928942509238843; 
    double xcrit = 192388.74826735;
    double havn = -1.0+1.0/(pow((T/xcrit),mu)+1.0);
    //Q = abs(havn);
    cout << "New Q = " << havn << endl;
//    cout << "pow(Q,2) = " << pow(Q,2) << endl;
//    cout << "q1 = " << setprecision(8) << q1 << endl;
//    cout << "q2 = " << q2 << endl;
    return (Q + (q1*Q*(Q-1.0)) + (q2*pow(Q,2)*(Q-1.0)));
  }


  SFunk::SFunk(const RDF& rdf, const SIM& sim, double T){
    
    rho = sim.get_rho();
    N = sim.get_N();
    cout << "N from constructor " << N << endl;
    mass = sim.get_mass();
    R = sim.get_R();
    cout << "mass from constructor " << mass << endl;
    r = rdf.get_r();
    gofr = rdf.get_gofr();
    ibins = rdf.get_ibins();
    cout << "ibins from constructor " << ibins << endl;

    //     De Broglie thermal wavelength calc and conversion
    dbw2 = pow(h,2.0)/(2.0*pi*mass*kB*T);
    dbw2 = dbw2*1e-15;
    dbw = sqrt(dbw2);
    dbw = dbw*1e10;
    
    //     1st Term - Sackur Tetrode Equation
    S_first = 2.5 - log(rho*pow(dbw,3.0));


    //  2nd Term - g(r)*ln[g(r)]*r^2 integral
    //  Trapezoidal Rule = 0.5*dx*(f(1) + 2*sum of f(2 - n-1) + f(n) )
    sum = 0.0;
    for (int i = 1; i<ibins-1; i++){
      if (gofr[i] > 0){
        sum = sum + gofr[i]*log(gofr[i])*pow(r[i],2);
      }
    }
    if (gofr[0] > 0){
      term_first = gofr[0]*log(gofr[0])*pow(r[0],2);
    }  
    if (gofr[ibins-1] > 0) {
      term_last= gofr[ibins-1]*log(gofr[ibins-1])*pow(r[ibins-1],2);
    }
    dr = r[1] - r[0];


    //  second S term
    S_second = 0.5*dr*(term_first + 2.0*sum + term_last);
    cout << "term_first " << term_first << endl;
    cout << "term_last " << term_last << endl;
    cout << "S_second " << S_second << endl;

    //     3rd Term - [g(r)-1]*r^2 integral
    //         Trapezoidal Rule = 0.5*dx*(f(1) + 2*sum of f(2 - n-1) + f(n) )
    sum = 0.0;
    
    for (int i = 1; i<ibins-1; i++){
      sum = sum + (gofr[i]-1.0)*pow(r[i],2);
    }

    term_first = (gofr[0]-1.0)*pow(r[0],2);
    term_last= (gofr[ibins-1]-1.0)*pow(r[ibins-1],2);

    S_third = 0.5*dr*(term_first + 2.0*sum + term_last);
    cout << "term_first " << term_first << endl;
    cout << "term_last " << term_last << endl;
    cout << "S_third " << S_third << endl;

    //      print *, 'third term', S_third
      

    //     Entropy (S_excess/N*kB, S/N*kB, and S)

    S_excess_over_NkB = 4.0*rho*pi*(S_second - S_third);
    S_over_NkB = S_first + S_excess_over_NkB;
    // S = S_over_NkB*N*kB;
    
    Q = Qofg(rho, N, r, gofr, ibins, dr, R); 
    cout << "Old Q = " << Q << endl;
    phi = phi_tilde(Q, q1, q2);
    cout << "phi_tilde = " << phi << endl;
    kappa_4 = pow(kappa(rho, r, gofr, ibins, dr),4);  
    cout << "1 + q0kappa_4 = " << setprecision(12) << (1.0 + q0*kappa_4) << endl;
    cout << "q0 " << q0 << endl;

  }

  //double SFunk::get_S() {
  //	  return S;
  //}

  double SFunk::get_S_excess_over_NkB() {
    return S_excess_over_NkB;
  }

  double SFunk::get_S_excess() {
 //   return -0.5*(1.0 + S_excess_over_NkB);
 //   return -0.5*(1.0 + S_excess_over_NkB)/(1.0 + (q0*kappa_4));
 //   return -0.5*(1.0 + S_excess_over_NkB)/(1.0 + (q0*kappa_4)) + 0.5*phi - 1.0;
   return -0.5*(-1.0 + S_excess_over_NkB) + 0.5*phi -1.0;
//    return (-0.5*(1.0 + S_excess_over_NkB))/(1.0 + (q0*kappa_4)) + 0.5*phi;

  }

  double SFunk::get_Q() {
    return Q;
  }

  
}
