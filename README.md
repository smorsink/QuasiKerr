# QuasiKerr with redshifts
Computes the RNS metric and the QuasiKerr Metric with various redshift values, with the option to compute redshift maps from different approximations for the zero-angular momentum-redshift, and using Oblate Schwarzschild light curves.

This is a set of C routines that compute the structure of a rapidly rotating star in order to compare its properties (including geodesics and redshift approximations) with a similar star computed in the quasikerr or oblate Schwarzschild approximations.

It also includes a Python file to compile, read, fit to and plot data sets computed with these routines, as well as the equations of state used in the data sets, in the 'compose' folder.

Redshifts are stored and plotted as 'z' instead of '1+z' or '1/(1+z)'

Usage:
1. Alter the makefile if you need to change the compiler. The computational gridsize is specified in the makefile. 
You can change to a smaller gridsize (like standard) to make the code run faster when you are just in the debugging
stage. For final results, I recommend the 201x401 grid. Compile code by typing at the command line:
> make
This creates an executable file named "quasik"

2. Choose an equation of state file, specify the file, the central energy density and the spin frequency or oblateness. 
I store the eos files in a directory called "compose" inside the directory with the code. So the path to 
my eos directory is "compose" You may have to change your path. As an example I'm running eos APR with a 
central energy density/c^2 of 2x10^{15} g/cm^3, and a spin frequency of 700 Hz. Type at the command line:
> quasik -f compose/eosAPR -e 1e15 -s 700
You may also specify the oblateness or "r ratio" instead of spin frequency using '-r'.

3. The code then spends about a minute finding the star with the correct parameters. The computation takes only a few seconds if the oblateness is specified instead. Once it finds the correct star, the rest of the computations (like of the metric for instance) are pretty fast. The output includes:
(a) Basic info about the star, like its mass, equatorial radius, and dimensionless angular momentum parameter.
(b) Info about the Quadrupole Moment using the Laarakkers and Poisson formula. The parameter epsilon is the same parameter as 
in the Glempadakis and Babak paper. epsilon = eta * a^2 where eta is the parameter in the paper by Baubock et al. Note that the Laarakkers and Poisson quadrupole moment is not exactly correct. The code also computes the correct value, as described by Pappas and Apostolatos. Look in the file "spin.c" for the definition of epsilon. You can comment out the Laarakkers and Poisson and switch to the other version if you like.
(c) The values of the radial coordinates at each value of lattitude on the star's surface are stored. 
This includes the isotropic, Schwarzschild, and Boyer-Lindquist values. You can use this to figure out appropriate values 
of radius to start your geodesic integrations at.
(d) The angular frequency and Kepler frequency

4. The code saves computed values in CSV files. To compute redshift maps as a function of angular coordinates, type '-m 1' when running the command line. (This increases computation time by a few seconds.) The maps are computed with a given inclination angle, specified in degrees using '-i', 90° by default. These files can be read and plotted using the Python notebook file 'plotting.ipynb'.
(a) The 'star.csv' file contains global parameters. These are the values saved in the data sets and used for the fits.
	[0] central density (1e15 g/cm^3), [1] angular velocity (rad/s), [2] Schwarzschild equatorial radius (cm), [3] mass (g), [4] M/R, [5] Omega_bar, [6] RNS equatorial b=0 redshift, [7] dimensionless quadrupole moment, [8] dimensionless angular momentum, [9] quadrupole correction factor epsilon, [10] fractional error of RNS and Oblate Scwarzschild redshift at the north pole, [11] Quasikerr equatorial b=0 redshift with fitted quadrupole moment
(b) The 'redshift.csv' file contains b=0 and b_extreme redshifts as a function of latitude.
	[0] mu, [1] RNS b_z=0, [2] RNS b_z=b_max, [3] RNS b_z=b_min, [4] Schwarzschild, [5] e^-(gamma+rho)/2, [6] Obl. Sch. b=0, [7] Obl. Sch. b=b_max, [8] Obl. Sch. b=b_min, [9] Lorentz gamma factor, [10] Kerr b=0, [11] Quasikerr b=0, [12] Quasikerr b=0 with fitted q
(c) The 'zmap.csv' file contains a 2D redshift map of the surface. The phi-dependence is computed the same way for each redshift, using the Oblate Schwarzschild approximation.
	[0] mu, [1] phi, [2] Obl. Sch, [3] RNS, [4] Quasikerr, [5] Quasikerr with fitted q

5. The files in the quasikerr directory are:

plotting.ipynb: Python file to compile and read data sets, compute fits, and make plots

star_data_EOS.npy: numpy files containing neutron star data sets for the specified EOS

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




