# SFunk
Code to support excess entropy calculations from classical molecular dynamics (MD) simulation trajectories following the methods presented in [Entropy Pair Functional Theory](https://doi.org/10.3390/e23020234)
The entropy functional has two forms, one for the liquid and gas phases (fluid) and one for the solid (crystal). Both are derived from the Kirkwood approximation for pair correlations. 


## Instalation
Clone this git repository
Navigate to either the solid or liquid source directory and edit the makefile to reflect your tooling setup
```
make
```

## Usage
A json configuration file is used to configure the software for your system
```
{
	"neighbors": 8,   // number of first nearest neighbors 
	"atoms": 93312,   // number of atoms per frame
	"frames": 10000,  // number of frames in the trajectory file 
	"skipframes": 0,  // number of frames to skip at the beginning of the trajectory
	"skin": 6.0,      // size in angrstroms of the padding to apply around the simulation cell to mimic periodic boundary conditions 
	"datafile": "Fe_dummy_long_2.trj", // path to the trajectory file to be read
	"outfile": "test.out",  // output file that will contain a record of the running covariance statistics (used to test for convergence of the Welford method)
	"dump": 10000,    // how many data points to skip when dumping running statistics into "outfile"
	"avgdump": false, // flag to set whether or not to write the average frame to a file (useful if your job runs out of time on an HPC resource)
	"avgoutfile": "avgout.trj", // file to write the "avgdump" to if "avgdump" is true
	"nbList": false,  // flag to set whether or not to write the neighbor list to a file (useful to examine local behavior of the system)
	"neighbor list file": "nbl.txt",  // file to write neighbor list to if "nbList" is true
	"avgtraj": "avg.trj", // hmmm, I forgot what this was for
	"var01_only": true,   // only perform covariance calculations taking the "avgoutfile" from a previous calculation 
	"avgtrajin": "avgout.trj" // a file containing a single average frame of data from a previous calculation
}
```
The code can then be run as
```
SFunk configfile.json > outputfilename
```





## Fluid Excess Entropy
The liquid/gas form of the functional takes as input a radial distribution function (RDF) that can be obtained from any MD software package or user developed post processing code. This form of the functional is not very computationally expensive and consequently a python package is provided in addition to the C++ code. At this time, for each new elemental system being studied three parameters will need to be fit in order to have a usuable functional for computing absolute entropy. A set of tools utilizing scipy curve fitting functions is provided as a starting point for the fitting of these parameters. Work is underway to eliminate the need to fit these parameters.


## Crystal Excess Entropy
The solid form of the functional takes as input a trajectory as output by MD simulations. The ordered nature of crytaline materials precludes the need for the calculation of an RDF. Instead autocovariance and first nearest neighbor covrainces are calculated directly from the positional data contained in the trajectory files. The [nanoflann](https://github.com/jlblancoc/nanoflann) header online kD-tree implementation is used to calculate a list of nearest neighbors based on each atom's ideal location in the system. For each atom in the system the ideal location is assumed to be the avarage position of the atom over the course of the simulation. This average position as well as the autocovariance is calculated in a single pass of the data. After the neighbor list has been generated the





@misc{blanco2014nanoflann,
  title        = {nanoflann: a {C}++ header-only fork of {FLANN}, a library for Nearest Neighbor ({NN}) with KD-trees},
  author       = {Blanco, Jose Luis and Rai, Pranjal Kumar},
  howpublished = {\url{https://github.com/jlblancoc/nanoflann}},
  year         = {2014}
}
