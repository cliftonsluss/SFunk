#include "sim.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace entropy_functional {

  SIM::SIM(double a, double box, int num_atoms, double atom_mass, double R_cut) {
  
    cout << "N from sim constructor " << N << endl; 
    V = pow(a*box,3);
    N = num_atoms;
    mass = atom_mass;
    rho = N/V;
    R = R_cut;
  }

  double SIM::get_rho() const{
    cout << "rho from get_rho() " << rho << endl;
    return rho;
  }

  double SIM::get_mass() const{
    cout << "mass from get_mass() " << mass << endl;
    return mass;
  }

  int SIM::get_N() const{
    cout << "N from get_N() " << N << endl;
    return N;
  }

  double SIM::get_R() const{
    cout << "R from get_R() " << R << endl;
    return R;
  }


}

