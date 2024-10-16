/************************************************************************** 
*                         quasikerr.c                                     *
*
* This program first computes the metric of a rotating star using rns.
* Then it computes the angular momentum parameter and quadrupole moment.
* It then computes three different analytic metrics:
*   (1) The Schwarzschild metric with the same value of M
*   (2) The Kerr metric with the same values of M and a
*   (3) The Quasi-Kerr metric with the same values of M, a, and q
*
* Makes use of the Quasi-Kerr metric defined in the paper:
* Glampedakis & Babak, "Mapping spacetimes with LISA: inspiral of a test
* body in a `quasi-Kerr' field", Classical and Quantum Gravity 23, 4167 (2006)
*
*                                                                         *
*                                                                         *
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "consts.h"
#include "struct.h"

#include "nrutil.h"
#include "equil.h"
#include "equil_util.h"
#include "findmodel.h"
#include "quadrupole.h"
#include "surface.h"
#include "interp.h"
#include "metric.h"


/* Main; where it all starts and ends */

int main(int argc, char **argv)     /* Number of command line arguments, Command line arguments */
{ NeutronStar star;
  EOS eos;
      
  int i, ierr;
  int m;

  double
    r_ratio=0.0,
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq, 
    Gamma_P=0.0,
    qscale
   ;

  
  double rbl, // Boyer-Lindquist radial coordinate
    gsch[5],  // Schwarzschild metric tensor components
    gmunu[5], // Kerr metric tensor components
    geps[5],  // Quasi-Kerr metric tensor components
    GammaSch[4][4][4], // Christoffel symbols for Schwarzschild
    Gamma[4][4][4],    // Christoffel symbols for Kerr
    GammaEps[4][4][4]; // Christoffel symbols for Quasi-Kerr
  double red_rns, red_sch0, red_sch, red_ker, red_quk; // Redshifts for rns, Schwarzschild, Kerr and Quasi-Kerr

  double mu, theta, epsilon;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  
  FILE *output;

  /* READ IN THE COMMAND LINE OPTIONS */
  for(i=1;i<argc;i++) 
    if(argv[i][0]=='-'){
      switch(argv[i][1]){

      case 'q':
	/* CHOOSE THE EOS TYPE: EITHER "tab" or "poly" or "quark"
	   (default is tab) */
	sscanf(argv[i+1],"%s",eos_type);
	break;  

      case 'b':
	sscanf(argv[i+1],"%lf",&B);
	B *= 1.602e33*KSCALE;
	break;       

      case 'f':
	/* IF A TABULATED EOS WAS CHOSEN, CHOOSE THE
	   NAME OF THE FILE */
	sscanf(argv[i+1],"%s",eos_file);
	break;

      case 'e':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_min);
	if(strcmp(eos_type,"poly")!=0)
	  e_min *= C*C*KSCALE;
	e_max = e_min;
	break;

      case 'l':
	/* CHOOSE THE CENTRAL ENERGY DENSITY OF THE 
	   NEUTRON STAR (IN g/cm^3) */
	sscanf(argv[i+1],"%lf",&e_max);
	if(strcmp(eos_type,"poly")!=0)
	  e_max *= C*C*KSCALE;
	break;

     case 's':
	/* CHOOSE THE SPIN FREQUENCY (HZ) */
	sscanf(argv[i+1],"%lf",&spin_freq);
	break;

      case 'r':
	/* CHOOSE r_ratio */
	sscanf(argv[i+1],"%lf",&r_ratio);
	break;

      case 'h': 
	fprintf(stderr,"\nQuick help:\n\n");
	fprintf(stderr,"  -q EOS type (tab)\n"); 
	fprintf(stderr,"     tab   : tabulated \n");
        fprintf(stderr,"     quark : simple quark model \n"); 
	fprintf(stderr,"  -b bag constant in MeV/fm^3 for quark models\n");
	fprintf(stderr,"  -f EOS file \n");
	fprintf(stderr,"  -e lowest central energy density to be used, in gr/cm^3\n");
	fprintf(stderr,"  -s spin frequency in Hz\n");
	fprintf(stderr,"  -r r_ratio [Default 0; if not 0 overrides spin frequency] \n");
	fprintf(stderr,"  -h this menu\n\n");
	exit(1);
	break;  
      }
    }


  /* PRINT THE HEADER */
  if(strcmp(eos_type,"tab")==0)
    printf("EOS file: %s, Grid size = MDIVxSDIV = %dx%d\n",eos_file,MDIV,SDIV);
  if(strcmp(eos_type,"quark")==0)
    printf("Quark star with B=%f, MDIVxSDIV=%dx%d\n",B/1.602e33/KSCALE,MDIV,SDIV);

  /* SetUpStar loads in the eos and sets up the grid */
  /* Source code for SetUpStar can be found in findmodel.c */

  ierr = SetUpStar(eos_file, eos_type, data_dir, Gamma_P, B, K,
		    &eos, &star);

  e_center = e_min;

    // Create a spherical star before computing the rotating star
    ierr = MakeSphere(&eos, &star, e_center);

 
    // Compute the structure of a rotating star
    if (r_ratio==0)
      ierr = SetSpin(&eos, &star, e_center, spin_freq);

    else
      ierr = rns(r_ratio, e_center, &eos, &star);

    qscale = quadrupole(&eos,&star);

    printf("\ne_center = %ge15 g/cm^3\n",star.e_center);
    printf("r_ratio  = %g\n", star.r_ratio);
    printf("Mass     = %g MSUN\n", star.Mass/MSUN);
    printf("M/R      = %g \n", star.Mass*G/(star.R_e*C*C));
    printf("R_e      = %g km\n",star.R_e*1e-5);

    //  printf("R_e      = %g [Mass units] \n", star.R_e*C*C/(star.Mass * G));
    //printf("R_iso    = %g km\n",star.r_e*sqrt(KAPPA)*1e-5);
    //printf("R_iso    = %g [Mass units] \n", star.r_e*sqrt(KAPPA)*C*C/(star.Mass*G));
    printf("Spin     = %g Hz\n",star.Omega/(2.0*PI));
    printf("j = a/M  = J/M^2 = %g \n\n", C/G * star.ang_mom/pow(star.Mass,2));
    //printf("j = a/R  = J/M^2 * M/R = %g \n\n", C/G * star.ang_mom/pow(star.Mass,2) * star.Mass*G/(star.R_e*C*C) );

    //printf("J = %g [g cm^2/s] \n",star.ang_mom);

    //printf("2J/R^3 = %g [g/cm 1/s] \n", 2.0*star.ang_mom * pow(star.R_e,-3));

    //printf("2J/R^3 = %g [code units] \n", 2.0*star.ang_mom * pow(star.R_e,-3) * G * sqrt(KAPPA)/(C*C*C));

    //printf("omega_surface = %g [code units] \n", star.metric.omega[(SDIV-1)/2+1][1]);
    //printf("Angular Veloc = %g rad/s \n", star.Omega);
    //printf("Angular Veloc = %g [code units] \n", star.Omega*sqrt(KAPPA)/C);


 
    printf("Quadrupole: Laarakkers and Poisson \n");
    printf("Quadrupole = %g g cm^2 \n", star.Quad);
    // qscale is the scaled quadrupole moment defined by Laarakkers and Poisson
    printf("qscale = %g \n",qscale);
    printf("mscale = %g \n",star.m_scale);
    //printf("qscale * (M/R)^3 = %g \n",qscale * pow(star.Mass*G/(star.R_e*C*C),3) );
    printf("LP formula: q = %g \n", - 4.3 * pow( C/G * star.ang_mom/pow(star.Mass,2) ,2));

    epsilon = star.epsilon;

    // Pappas and Apostolatos definition
    //epsilon = -star.m_scale - pow(star.j,2);
    printf("QuasiKerr dimensionless quadrupole parameter: epsilon = %g \n",epsilon);

    
    // Find the Star's surface. 
    ierr = Surface(&eos,&star);
    
    //SurfPrint(&eos,&star); Uncomment this line if you'd like information about the surface printed out to a file.

    // Set up the Boyer-Linquist Coordinate system.
    // Code can be found in quadrupole.c
    ierr = BLSetUp(&star);

    //printf("R_S      = %g [Mass units] \n", star.blmetric.r_S[(SDIV-1)/2+1]);
    
    // Example of how to make make the function calls and use of the metric and Christoffel symbols
    // Choose some values of r_BL and theta and then evaluate!
    printf("\nComparison of QuasiKerr Metric with numerical RNS metric\n");

    output = fopen("redshift.txt","w");
    fprintf(output, "#Redshift values on surface (Energy Obs/Energy Emitted) \n"
	    );
    fprintf(output, "#theta      mu         RNS   Schw(0)   Schw   Kerr    Quasi-K   \n"
	    );

    for (m=1;m<MDIV;m++){
      
      mu = star.metric.mu[m];
      theta = acos(mu);

      // Choose the value of Boyer-Linquist radial coordinate on the surface. 
      rbl =  star.r_BL_surf[m];   

      // This function call computes transforms the RNS metric into the BL coordinates at this point
      ierr = BLComputeMetric(&star, rbl, theta);

      //if (m==1)
      //	printf(" (v/c) = %g   (v/c)^2 = %g = %g \n", star.Omega/C * star.R_e,  pow(star.Omega/C * star.R_e,2),
      //	  pow(star.Omega*sqrt(KAPPA)/C,2) * star.blmetric.g_BL[3] * pow(star.Mass * G/(C*C*sqrt(KAPPA)) ,2)   );

      // Michi's routines for computing the QuasiKerr metric

      // metric computes the value of the quasi-kerr metric at a given value of rbl and theta
      // function requires value of the Kerr parameter a = star.j and the epsilon parameter
      
      // gsch stores the metric with a = epsilon = 0 -- Schwarzschild
      metric(rbl, theta, gsch, 0, 0);

      // gmunu stores the metric with epsilon = 0 -- Kerr
      metric(rbl, theta, gmunu, star.j, 0);

      // geps stores the metric with nonzero epsilon
      metric(rbl, theta, geps,star.j,epsilon);

      

      /*	if (m==1){
	  printf("Comparison of Numerical RNS metric with QuasiKerr metric approx. \n");
	  printf("         RNS  Schw    Kerr   QK \n");
	   printf("gtt   = %g   = %g    = %g   = %g \n",
	   star.blmetric.g_BL[0],gsch[0],gmunu[0],geps[0]);
	   printf("grr   = %g   = %g    = %g   = %g \n",
		  star.blmetric.g_BL[1],gsch[1],gmunu[1],geps[1]);
	   printf("gthth = %g   = %g    = %g   = %g \n",
		  star.blmetric.g_BL[2],gsch[2], gmunu[2],geps[2]);
	   printf("gpiph = %g   = %g    = %g   = %g \n",
		  star.blmetric.g_BL[3],gsch[3], gmunu[3],geps[3]);
	   printf("gphit = %g   = %g    =%g    =%g\n",
		  star.blmetric.g_BL[4],gsch[4], gmunu[4],geps[4]);
	   printf("gtphi/gphiphi = %g \n",
		  star.blmetric.g_BL[4]/star.blmetric.g_BL[3]/(star.Mass * G/(C*C*sqrt(KAPPA))));
		  }*/

      // star.blmetric.g_BL = metric function for RNS
      // gmunu = Kerr
      // geps = quasi-kerr
      // Definitions of metric components:
      // gmunu[0]=g_tt
      // gmunu[1]=g_rr
      // gmunu[2]=g_thth
      // gmunu[3]=g_phiphi
      // gmunu[4]=g_phit
      

	// The metric components gmunu[3] are written in dimensionless "mass units"
	// The Angular velocity must be converted to dimensionless mass units. To do this, note that angular velocity has
	// units of 1/time. Meanwhile GM/c^3 has units of time. Therefore to convert, multiply every instance
	// of Omega by GM/c^3.

      // redshift factors are the ratio E_obs/E_em
      // E_obs = photon energy measured by an observer at infinity
      // E_em = photon energy emitted at the star's surface
      
	// red_sch0 is the redshift factor WITHOUT the transverse doppler effect
      red_sch0 = pow( -1.0 * (gsch[0] 
			     ), 0.5);

      // red_sch includes the transverse doppler effect
      red_sch = pow( -1.0 * (gsch[0] 
			     + pow( star.Omega * G/(C*C*C) * star.Mass,2 ) * gsch[3] 
			     ), 0.5);
      

      red_rns = pow( -1.0 * (star.blmetric.g_BL[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * star.blmetric.g_BL[4]
			     + pow( star.Omega * G/(C*C*C) * star.Mass,2 ) * star.blmetric.g_BL[3]
			     ), 0.5);

      red_ker = pow( -1.0 * (gmunu[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * gmunu[4] 
			     + pow( star.Omega * G/(C*C*C) * star.Mass,2 ) * gmunu[3]
			     ), 0.5);

      red_quk = pow( -1.0 * (geps[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * geps[4]
			     + pow( star.Omega * G/(C*C*C) * star.Mass,2 ) * geps[3] 
			     ), 0.5);
      
      
      if (m==1){
	printf("Sch (1+z) = %g  RNS (1+z) = %g \n", red_sch, red_rns);
	printf(" Fractional diff in Flux = %g\n", pow(red_sch/red_rns,3)-1.0);
      }

      fprintf(output, "%d %lf %lf  %lf %lf %lf %lf %lf \n",
	      m,theta, mu, red_rns, red_sch0, red_sch, red_ker,red_quk);
     

    }





    
    /*

    // Michi's Christoffel Symbol routines
    //Gamma stores the Christoffel symbols for epsilon = 0
    gammas(rbl, theta, Gamma, star.j, 0.0);
    // Gamma store the Christoffel symbols for nonzero epsilon
    gammas(rbl, theta, GammaEps, star.j, epsilon);
    // GammaSch store the Christoffel symbols for Schwarzschild
    gammas(rbl, theta, GammaSch, 0, 0);

    printf("Comparison of Christoffel Symbols\n");

    for(i=0;i<=3;i++)
      for(j=0;j<=3;j++)
	for(k=0;k<=3;k++)
	  printf("Gamma[%d][%d][%d] = %g = %g = %g %g \n",i,j,k,
		 star.blmetric.Gamma_BL[i][j][k], Gamma[i][j][k], GammaEps[i][j][k], GammaSch[i][j][k]);


    */

  return 0;
}









