#include "sfunk_const.h"

#ifndef SIM_H
#define SIM_H

namespace entropy_functional {

  class SIM {

    public:
      SIM(double a, double box, int num_atoms, double atom_mass, double R_cut);
      double get_rho() const;
      double get_mass() const;
      int get_N() const;
      double get_R() const;
      void setup(int N, double mass);

    private:
      double V;
      double rho;
      double mass;
      double R;
      int N;

  }; 

}
#endif


