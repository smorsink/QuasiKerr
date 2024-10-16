# QuasiKerr
Computes the RNS metric and the QuasiKerr Metric

This is a set of C routines that compute the structure of a rapidly rotating star in order to compare its properties (including geodesics) with a similar star computed in the quasikerr approximation. 

Usage:
1. Alter the makefile if you need to change the compiler. The computational gridsize is specified in the makefile. 
You can change to a smaller gridsize (like standard) to make the code run faster when you are just in the debugging
stage. For final results, I recommend the 201x401 grid. Compile code by typing at the command line:
> make
This creates an executable file named "quasik"

2. Choose an equation of state file, specify the file, the central energy density and the spin frequency. 
I store the eos files in a directory called "eos" parallel to the directory with the code. So the path to 
my eos directory is "../eos" You may have to change your path. As an example I'm running eos APR with a 
central energy density/c^2 of 2x10^{15} g/cm^3, and a spin frequency of 700 Hz. Type at the command line:
> quasik -f ../eos/eosAPR -e 2e15 -s 700

A simple bash script go.sh is included. This script compiles the code an runs the code for EOS FPS, creating a star that spins at 700 Hz wiht a mass of 1.4 Msun. This is to help compare with the results of the Baubock et al 2013 paper. To run the script:
> go.sh

3. The code then spends about a minute finding the star with the correct parameters. Once it finds the correct star, the rest of the computations (like of the metric for instance) are pretty fast. The output includes:
(a) Basic info about the star, like its mass, equatorial radius, and dimensionless angular momentum parameter.
(b) Info about the Quadrupole Moment using the Laarakkers and Poisson formula. The parameter epsilon is the same parameter as 
in the Glempadakis and Babak paper. epsilon = eta * a^2 where eta is the parameter in the paper by Baubock et al. Note that the Laarakkers and Poisson quadrupole moment is not exactly correct. The code also computes the correct value, as described by Pappas and Apostolatos. Look in the file "spin.c" for the definition of epsilon. You can comment out the Laarakkers and Poisson and switch to the other version if you like.
(c) The values of the radial coordinates at each value of lattitude on the star's surface are stored. 
This includes the isotropic, Schwarzschild, and Boyer-Lindquist values. You can use this to figure out appropriate values 
of radius to start your geodesic integrations at.
(d) The program then chooses a value of radius and co-lattitude and computes the value of the metric and Christoffel symbols. The redshift for zero angular momentum photons is computed using the various metrics and output to a file named redshift.txt.

4. The files in the quasikerr directory are:

quasikerr.c = main routine

consts.h = constants used by the various routines

equil.c, equil.h, equil_util.c, equil_util.h = routines used to compute the rapidly rotating neutron star

findmodel.c, findmodel.h = routines used to control the star's parameters (like the spin)

interp.c, interp.h = various interpolation routines

makefile = compilation instructions

metric.c, metric.h = the QuasiKerr metric code that Michi Baubock wrote

nrutil.c, nrutil.h = various routines from Numerical Recipes

quadrupole.c, quadrupole.h = routines for computing the quadrupole moment as well as the RNS metric and
	      Christoffel symbols at any point outside of the star.
       
struct.h = definitions for all the structures used in the code

surface.c, surface.h = routines that find the star's surface

References

Laarakkers and Poisson,  "Quadrupole Moments of Rotating Neutron Stars", ApJ 512, 282, 1999

Glampedakis and Babak, "Mapping spacetimes with LISA: inspiral of a test body in a `quasi-Kerr' field", Class. Quatum Grav. 23, 4167, 2006

Baubock, Psaltis, and Ozel, "Narrow Atomic Features from Rapidly Spinning Neutron Stars", ApJ 766, 87, 2013

Pappas and Apostolatos, "Revising the Multipole Moments of Numerical Spacetimes and its Consequences", PRL 108, 231104, 2012




