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
