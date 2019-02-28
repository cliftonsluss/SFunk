#include "SFunk.h"
#include "rdf.h"
#include "sim.h"
#include <vector>
#include <iostream>

using namespace std;
using namespace entropy_functional;

int main() {


  	//double T = 1000000.0;
	double a = 2.86997;
	double T = 1000000.0;
	double box = 20.0;
	int N = 16000;
	//double R = 28.5995;
	double R = 28.5995;
        //double R = 28.60;
	double mass = 3.36299e-21; 
	string filename = "./rdfs/iron_1000000K.rdf";
	RDF rdf(filename);
	rdf.test();
	SIM sim(a, box, N, mass, R);
	SFunk s_funk(rdf, sim, T);

	cout << "S excess  = " << s_funk.get_S_excess() << "\n";
        cout << "Q = " << s_funk.get_Q() << "\n";
	cout << "rdf integral = " << rdf.get_integral() << endl;

}
