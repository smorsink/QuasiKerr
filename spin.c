/************************************************************************** 
*                         spin.c                                       *
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
      
  int i,j,k, ierr;
  int m,s;

  double
    r_ratio=0.0,
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq, 
    Gamma_P,
    qscale
   ;
   
  double rbl, riso, rs, gsch[5], gmunu[5], geps[5], GammaSch[4][4][4], Gamma[4][4][4], GammaEps[4][4][4];
  double red_rns, red_sch, red_ker, red_quk;

  double mu, theta, epsilon;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  
  FILE *output, *output2;

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

    printf("R_e      = %g [Mass units] \n", star.R_e*C*C/(star.Mass * G));
    printf("R_iso    = %g km\n",star.r_e*sqrt(KAPPA)*1e-5);
    printf("R_iso    = %g [Mass units] \n", star.r_e*sqrt(KAPPA)*C*C/(star.Mass*G));
    printf("Spin     = %g Hz\n",star.Omega/(2.0*PI));
    printf("j = a/M  = J/M^2 = %g \n\n", C/G * star.ang_mom/pow(star.Mass,2));
    printf("j = a/R  = J/M^2 * M/R = %g \n\n", C/G * star.ang_mom/pow(star.Mass,2) * star.Mass*G/(star.R_e*C*C) );

    printf("J = %g [g cm^2/s] \n",star.ang_mom);

    printf("2J/R^3 = %g [g/cm 1/s] \n", 2.0*star.ang_mom * pow(star.R_e,-3));

    printf("2J/R^3 = %g [code units] \n", 2.0*star.ang_mom * pow(star.R_e,-3) * G * sqrt(KAPPA)/(C*C*C));

    printf("omega_surface = %g [code units] \n", star.metric.omega[(SDIV-1)/2+1][1]);
    
 
    printf("Quadrupole: Laarakkers and Poisson \n");
    printf("Quadrupole = %g g cm^2 \n", star.Quad);
    printf("qscale = %g \n",qscale);
    printf("qscale * (M/R)^3 = %g \n",qscale * pow(star.Mass*G/(star.R_e*C*C),3) );
    
    // Find the Star's surface. 
    ierr = Surface(&eos,&star);
    SurfPrint(&eos,&star);

    // Set up the Boyer-Linquist Coordinate system.
    ierr = BLSetUp(&star);

    printf("R_S      = %g [Mass units] \n", star.blmetric.r_S[(SDIV-1)/2+1]);
    
    // Example of how to make make the function calls and use of the metric and Christoffel symbols
    // Choose some values of r_BL and theta and then evaluate!
    printf("Comparison of QuasiKerr Metric with numerical RNS metric\n");


    output = fopen("framedrag_mu.txt","w");
    //output = fopen("redshift.txt","w");
    fprintf(output, "#g_phit values on surface \n"	    );
    fprintf(output, "#theta   mu    RNS    Kerr    Quasi-K    (RNS-K)/RNS\n"	    );

    for (m=1;m<MDIV;m++){
      // Choose some value of theta
      //m = MDIV-1;
      mu = star.metric.mu[m];
      theta = acos(mu);
      // printf("theta = %g radians; mu = %g \n", theta, mu);

      // Choose the value of Boyer-Linquist radial coordinate on the surface. 
      rbl =  star.r_BL_surf[m];   
      //   printf(" s = %g \n", star.s_surf[m]);
      //printf("Boyer-Lindquist radial coordinate: rbl = %g M\n",rbl);

      // Laarakkers and Poisson definition
      epsilon = star.epsilon;
      // Pappas and Apostolatos definition
      //epsilon = -star.m_scale - pow(star.j,2);
      // printf("QuasiKerr dimensionless quadrupole parameter: epsilon = %g \n",epsilon);

      // This function call computes the RNS numerical metric and Christoffel symbols 
      ierr = BLComputeMetric(&star, rbl, theta);

      // Michi's routines for computing the QuasiKerr metric

      // gsch stores the metric with a = epsilon = 0 -- Schwarzschild
      metric(rbl, theta, gsch, 0, 0);

      // gmunu stores the metric with epsilon = 0 -- Kerr
      metric(rbl, theta, gmunu, star.j, 0);

      // geps stores the metric with nonzero epsilon
      metric(rbl, theta, geps,star.j,epsilon);
      /*
	printf("Comparison of Numerical RNS metric with QuasiKerr metric approx. \n");

	printf("gtt = %g = %g = %g = %g \n",
	star.blmetric.g_BL[0],gsch[0],gmunu[0],geps[0]);
	printf("grr = %g = %g = %g \n",
	star.blmetric.g_BL[1],gmunu[1],geps[1]);
	printf("gthth = %g = %g = %g \n",
	star.blmetric.g_BL[2],gmunu[2],geps[2]);
	printf("gphiphi = %g = %g = %g\n",
	star.blmetric.g_BL[3],gmunu[3],geps[3]);
	printf("gphit = %g = %g =%g \n",
	star.blmetric.g_BL[4],gmunu[4],geps[4]);
      */

      // star.blmetric.g_BL = metric function for RNS
      // gmunu = Kerr
      // geps = quasi-kerr
    
      fprintf(output, "%lf %lf %lf %lf %lf %lf \n",
	      theta, mu, star.blmetric.g_BL[4], gmunu[4], geps[4],(star.blmetric.g_BL[4]-gmunu[4])/star.blmetric.g_BL[4]);
    
    
      red_rns = pow( -1.0 * star.blmetric.g_BL[3]/(star.blmetric.g_BL[0]*star.blmetric.g_BL[3]-pow(star.blmetric.g_BL[4],2)),0.5);
      red_sch = pow( -1.0 * gsch[0],-0.5);
      red_ker = pow( -1.0 * gmunu[3]/(gmunu[0]*gmunu[3]-pow(gmunu[4],2)),0.5);
      red_quk = pow( -1.0 * geps[3]/(geps[0]*geps[3]-pow(geps[4],2)),0.5);
      /*
	printf("Static Redshift 1+z = (-g_{phi,phi}/(gttxg_{phi,phi}-g_{t,phi}^2))^{1/2}\n");
	printf("RNS 1+z = %lf \n",
	red_rns);

	printf("Schwarzschild 1+z = %lf, percent dif = %lf \n",
	red_sch, (red_rns-red_sch)/red_rns * 100.0 );

	printf("Kerr 1+z = %lf percent dif = %lf\n",
	red_ker, (red_rns-red_ker)/red_rns * 100.0);

	printf("Quasi-Kerr 1+z = %lf percent dif = %lf\n",
	red_quk, (red_rns-red_quk)/red_rns * 100.0);


	printf("Ang Vel = %lf  Frame-dragging = %lf \n",
	   star.Omega * sqrt(KAPPA)/C, star.blmetric.g_BL[4]/star.blmetric.g_BL[3]);

	   printf("Trans Dopp RNS = (1-v^2)^{-1/2} = %lf \n",
	   1.0/sqrt( 1.0 -  pow(star.blmetric.g_BL[3]*star.Omega * sqrt(KAPPA)/C 
	   + star.blmetric.g_BL[4],2.0)/fabs(star.blmetric.g_BL[0]*star.blmetric.g_BL[3]-pow(star.blmetric.g_BL[4],2))));
      */

      /*    printf("Trans Dopp Sch = (1-v^2)^{-1/2} = %lf \n",
	    1.0/sqrt( 1.0 -  pow(gsch[3]*star.Omega * sqrt(KAPPA)/C,2.0)/fabs(gsch[0]*gsch[3])));

	    printf("Trans Dopp Kerr = (1-v^2)^{-1/2} = %lf \n",
	    1.0/sqrt( 1.0 -  pow(gmunu[3]*star.Omega * sqrt(KAPPA)/C + gmunu[4],2.0)/fabs(gmunu[0]*gmunu[3]-pow(gmunu[4],2))));


	    printf("Trans Dopp QK = (1-v^2)^{-1/2} = %lf \n",
	    1.0/sqrt( 1.0 -  pow(geps[3]*star.Omega * sqrt(KAPPA)/C + geps[4],2.0)/fabs(geps[0]*geps[3]-pow(geps[4],2))));
      */


      //  fprintf(output, "%lf %lf %lf %lf %lf %lf \n",
      //	    theta, mu, red_rns, red_sch, red_ker, red_quk);
    }


   // Example of how to make make the function calls and use of the metric and Christoffel symbols
    // Choose some values of r_BL and theta and then evaluate!
    printf("Comparison of QuasiKerr Metric with numerical RNS metric\n");


    output = fopen("framedrag_r.txt","w");
    //output = fopen("redshift.txt","w");
    fprintf(output, "#g_phit values on equator \n"
	    );

    fprintf(output, "#r RNS Kerr Quasi-K diff\n"
	    );

    output2 = fopen("metric-equator.txt", "w");
    fprintf(output2, "#Values of Metric Potentials in the Equatorial Plane\n");
    fprintf(output2, "#Values of riso and r given in unit of M\n");
    fprintf(output2, "#Metric potentials rho, gamma, alpha are unitless \n");
    fprintf(output2, "#Metric potential has units of 1/length in rns code units. \n");
    fprintf(output2, "#riso    r        omega     2J/r^3    rho       gamma       alpha    -2M/r  \n");

    for (s=(SDIV-1)/2+1;s<SDIV;s++){


    // Choose some value of theta
      m=1;
      mu = star.metric.mu[m];
      theta = acos(mu);
      // printf("theta = %g radians; mu = %g \n", theta, mu);

    // Choose the value of Boyer-Linquist radial coordinate
      rbl =  star.blmetric.r_BL[s];
      riso = star.blmetric.r_iso[s];
      rs = star.blmetric.r_S[s];
      // printf(" s = %g \n", star.s_surf[m]);
    //printf("Boyer-Lindquist radial coordinate: rbl = %g M\n",rbl);

    /*  fprintf(output2,"%lf %lf %lf %g %g %g \n",
	      riso, rs, rbl, star.metric.omega[s][m],
	      2.0*star.ang_mom * pow(star.R_e,-3) * G * sqrt(KAPPA)/(C*C*C) *
	      pow( star.blmetric.r_BL[(SDIV-1)/2+1]/rbl   ,3),
	      2.0*star.ang_mom * pow(star.R_e,-3) * G * sqrt(KAPPA)/(C*C*C) *
	      pow( star.blmetric.r_S[(SDIV-1)/2+1]/rs   ,3));*/


 fprintf(output2,"%lf %lf %g %g %g %g %g %g \n",
	 riso, rs, star.metric.omega[s][m],
	 2.0*star.ang_mom * pow(star.R_e,-3) * G * sqrt(KAPPA)/(C*C*C) *
	 pow( star.blmetric.r_S[(SDIV-1)/2+1]/rs   ,3),
	 star.metric.rho[s][m],
	 star.metric.gama[s][m],
	 star.metric.alpha[s][m],
	 -2.0/rs
	 );
      

      
    // Laarakkers and Poisson definition
    epsilon = star.epsilon;
    // Pappas and Apostolatos definition
    //epsilon = -star.m_scale - pow(star.j,2);

    
    // printf("QuasiKerr dimensionless quadrupole parameter: epsilon = %g \n",epsilon);

    // This function call computes the RNS numerical metric and Christoffel symbols 
    ierr = BLComputeMetric(&star, rbl, theta);

    // Michi's routines for computing the QuasiKerr metric

    // gsch stores the metric with a = epsilon = 0 -- Schwarzschild
    metric(rbl, theta, gsch, 0, 0);

    //printf("star.j = %lf \n", star.j);
    // gmunu stores the metric with epsilon = 0 -- Kerr
    metric(rbl, theta, gmunu, star.j, 0);

    // geps stores the metric with nonzero epsilon
    metric(rbl, theta, geps,star.j,epsilon);

    //metric(rbl, theta, geps,0.08664,0);
    
    /*
    printf("Comparison of Numerical RNS metric with QuasiKerr metric approx. \n");

    printf("gtt = %g = %g = %g = %g \n",
	   star.blmetric.g_BL[0],gsch[0],gmunu[0],geps[0]);
    printf("grr = %g = %g = %g \n",
	   star.blmetric.g_BL[1],gmunu[1],geps[1]);
    printf("gthth = %g = %g = %g \n",
	   star.blmetric.g_BL[2],gmunu[2],geps[2]);
    printf("gphiphi = %g = %g = %g\n",
	   star.blmetric.g_BL[3],gmunu[3],geps[3]);
    printf("gphit = %g = %g =%g \n",
	   star.blmetric.g_BL[4],gmunu[4],geps[4]);
    */

    // g_phi_t
      fprintf(output, "%lf %lf %lf %lf %lf %lf %lf %lf %lf \n",
	   rbl, star.blmetric.g_BL[4], gmunu[4], geps[4],(star.blmetric.g_BL[4]-gmunu[4])/star.blmetric.g_BL[4],
	      (star.blmetric.g_BL[4]-geps[4])/star.blmetric.g_BL[4], riso, gmunu[4] * pow(riso/rbl,2),
	      (star.blmetric.g_BL[4] - gmunu[4] * pow(riso/rbl,2))/(star.blmetric.g_BL[4]) );

   // g_theta_theta
   /*  fprintf(output, "%lf %lf %lf %lf %lf %lf \n",
	   rbl, star.blmetric.g_BL[2], gmunu[2], geps[2],(star.blmetric.g_BL[2]-gmunu[2])/star.blmetric.g_BL[2],
	   (star.blmetric.g_BL[2]-geps[2])/star.blmetric.g_BL[2]);*/
 
       // g_phi_phi
    /*   fprintf(output, "%lf %lf %lf %lf %lf %lf \n",
	   rbl, star.blmetric.g_BL[3], gmunu[3], geps[3],(star.blmetric.g_BL[3]-gmunu[3])/star.blmetric.g_BL[3],
	   (star.blmetric.g_BL[3]-geps[3])/star.blmetric.g_BL[3]);*/

     // g_tt
    /* fprintf(output, "%lf %lf %lf %lf %lf %lf \n",
	   rbl, star.blmetric.g_BL[0], gmunu[0], geps[0],(star.blmetric.g_BL[0]-gmunu[0])/star.blmetric.g_BL[0],
	   (star.blmetric.g_BL[0]-geps[0])/star.blmetric.g_BL[0]);*/


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









