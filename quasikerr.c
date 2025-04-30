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
  int p;

  double
    r_ratio=0.0,
    e_min, e_max,
   e_center=1e15,                     /* central en. density */
   B,                            /* Quark Bag Constant */
   K=3.0,                        /* Second parameter in "quark" eos */
    spin_freq, 
    Gamma_P=0.0,
    incl_deg=90,                       /* observer inclination angle */
    eta_bad=3.3,                           /* wrong quadrupole correction */
    make_map=0,               /* 1 to make redshift map, 0 to skip */
    qscale
   ;

  
  double rbl, // Boyer-Lindquist radial coordinate
    gsch[5],  // Schwarzschild metric tensor components
    gmunu[5], // Kerr metric tensor components
    geps[5],  // Quasi-Kerr metric tensor components
    geps_fit[5],
    GammaSch[4][4][4], // Christoffel symbols for Schwarzschild
    Gamma[4][4][4],    // Christoffel symbols for Kerr
    GammaEps[4][4][4]; // Christoffel symbols for Quasi-Kerr
  double red_rns, red_sch0, red_sch, red_ker, red_quk, red_quk_fit; // Redshifts for rns, Schwarzschild, Kerr and Quasi-Kerr

  double mu, theta, epsilon;

  char eos_file[80] = "no EOS file specified";   /* EOS file name */
  char eos_type[80] = "tab";                     /* EOS type (poly or tab) */
  char data_dir[80] = "junk";                    /* Data output directory */
  //char map_buffer[1000];
  //size_t bufsize = 8*sizeof(double);
  
  FILE *output;

  // Stores b=0 and limb redshifts as a function of mu
  FILE *z_out;
  z_out = fopen("redshift.csv", "w");

  // Stores global parameters
  FILE *star_output;
  star_output = fopen("star.csv", "w");

  // Stores redshift maps
  FILE *z_map;
  z_map = fopen("zmap.csv", "w");

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

          case 'i':
	/* CHOOSE INCLINATION */
	sscanf(argv[i+1],"%lf",&incl_deg);
	break;

          case 'p':
  /* CHOOSE QUADRUPOLE CORRECTION */
  // currently not used
  sscanf(argv[i+1],"%lf",&eta_bad);
  break;

          case 'm':
  /* CHOOSE whether to compute redshift map */
  sscanf(argv[i+1],"%lf",&make_map);
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
  fprintf(stderr,"  -i inclination degrees (90)\n");
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
    //printf("\nComparison of QuasiKerr Metric with numerical RNS metric\n");

    /*
    output = fopen("redshift.txt","w");
    fprintf(output, "#Redshift values on surface (Energy Obs/Energy Emitted) \n"
	    );
    fprintf(output, "#theta      mu         RNS   Schw(0)   Schw   Kerr    Quasi-K   \n"
	    );*/

    printf("Omega = %lf, Omega_K = %lf\n", star.Omega, star.Omega_K);

    double z_00 = redshift_metric(&star, 0, 1);

    // find the approximated quadrupole moment and use it for QK b=0 redshift
    double quad_f = quad_fit(&star);
    double eps_f = -quad_f-SQ(star.j);
    //printf("q fit = %lf, q = %lf, eps fit = %lf, eps = %lf\n", quad_f, star.q, eps_f, star.epsilon);
    rbl = star.r_BL_surf[1];
    mu = star.metric.mu[1];
    theta = acos(mu);
    metric(rbl, theta, geps, star.j, eps_f);
    double z_qkf0 = pow( -1.0 * (geps[0] +
      2 * star.Omega * G/(C*C*C) * star.Mass * geps[4]
      + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * geps[3] 
      ), 0.5);

    //printf("RNS redshift = %lf, QK fit redshift = %lf\n", z_00, (1/z_qkf0)-1);

    double z_pole_error = (z_00-redshift_OS(&star, incl_deg, 0, 0, incl_deg*PI/180, MDIV))/(1+z_00);
    fprintf(star_output, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", 
    star.e_center*1e15, star.Omega, r_surf_sch(&star, 1), star.Mass, star.Mass*G/(SQ(C)*r_surf_sch(&star, 1)), star.Omega*sqrt(r_surf_sch(&star, 1)*SQ(r_surf_sch(&star, 1))/(star.Mass*G)), z_00, star.q, star.j, star.epsilon, z_pole_error, (1/z_qkf0)-1);
    fclose(star_output);

    if (make_map == (double)1){
      printf("Making redshift map\n");
    double incl = incl_deg*PI/180;
    double *b = (double *)malloc((MDIV+1) * sizeof(double));
    double *psi_b = (double *)malloc((MDIV+1) * sizeof(double));
    double *zqk_arr = (double *)malloc((MDIV+1) * sizeof(double));
    double *zqkb_arr = (double *)malloc((MDIV+1) * sizeof(double));
    double nothing = 0;
    if (b == NULL || psi_b == NULL || zqk_arr == NULL || zqkb_arr == NULL) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(EXIT_FAILURE);
    }
    for (m=MDIV;m>1;m--){
      
      mu = star.metric.mu[m];
      theta = acos(mu);

      // Choose the value of Boyer-Linquist radial coordinate on the surface. 
      rbl =  star.r_BL_surf[m];   

      // This function call computes transforms the RNS metric into the BL coordinates at this point
      ierr = BLComputeMetric(&star, rbl, theta);
      
      // gsch stores the metric with a = epsilon = 0 -- Schwarzschild
      metric(rbl, theta, gsch, 0, 0);

      // gmunu stores the metric with epsilon = 0 -- Kerr
      metric(rbl, theta, gmunu, star.j, 0);

      // geps stores the metric with nonzero epsilon
      metric(rbl, theta, geps, star.j, epsilon);

      metric(rbl, theta, geps_fit, star.j, -quad_fit(&star)-SQ(star.j));

      //fprintf(surface_out, "%lf, %lf \n", star.metric.mu[m], r_surf_sch(&star, m)*1e-5);
      
	// red_sch0 is the redshift factor WITHOUT the transverse doppler effect
      red_sch0 = pow( -1.0 * (gsch[0] 
			     ), 0.5);

      // red_sch includes the transverse doppler effect
      red_sch = pow( -1.0 * (gsch[0] 
			     + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * gsch[3] 
			     ), 0.5);
      

      red_rns = pow( -1.0 * (star.blmetric.g_BL[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * star.blmetric.g_BL[4]
			     + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * star.blmetric.g_BL[3]
			     ), 0.5);

      red_ker = pow( -1.0 * (gmunu[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * gmunu[4] 
			     + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * gmunu[3]
			     ), 0.5);

      red_quk = pow( -1.0 * (geps[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * geps[4]
			     + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * geps[3] 
			     ), 0.5);

      red_quk_fit = pow( -1.0 * (geps_fit[0] +
			     2 * star.Omega * G/(C*C*C) * star.Mass * geps_fit[4]
			     + pow(star.Omega * G/(C*C*C) * star.Mass, 2) * geps_fit[3] 
			     ), 0.5);   

      zqk_arr[m] = red_quk;
      zqkb_arr[m] = red_quk_fit;
      
      
      if (m==1){
	      printf("Sch (1+z) = %g  RNS (1+z) = %g \n", red_sch, red_rns);
	      printf(" Fractional diff in Flux = %g\n", pow(red_sch/red_rns,3)-1.0);
      }

      double z_0 = redshift_metric(&star, 0, m);
      double z_max = redshift_metric(&star, b_extreme_metric(&star, (double)1, m), m);
      double z_min = redshift_metric(&star, b_extreme_metric(&star, (double)-1, m), m);

      double z_sch = redshift_schwarzschild(&star, m);
      double z_grav = redshift_grav_metric(&star, m);
      
      double b_max_OS = b_extreme_OS(&star, m);
      
      double z_os_0 = redshift_OS(&star, incl_deg, 0, 0, 0, m);

      double gam = gamma(&star, m);

      // compute map
      // northern hemisphere
      
      double b_step = (2*b_max_OS)/(MDIV-2);
      double phi_step = (2*PI)/(MDIV-2);
      double psi_b_max = psi_integrated(&star, b_max_OS, m);
      double psi_b_min = psi_integrated(&star, -b_max_OS, m);
      double b_final_max=0, b_final_min=0, phi_max, phi_min, psi_max, psi_min;
      for (p=1;p<MDIV;p++){
        
        double b_trial = -b_max_OS+(p-1)*b_step;
        double psi_trial = psi_integrated(&star, b_trial, m);

        b[p] = b_trial;
        psi_b[p] = psi_trial;
      }
      //fprintf(z_map, "%lf\n", star.metric.mu[m]);
      for (p=1;p<MDIV;p++){

        double phi = -PI+(p-1)*phi_step;

        double psi = acos(cos(incl)*star.metric.mu[m]+sin(incl)*sqrt(1-SQ(star.metric.mu[m]))*cos(phi));

        if (psi > psi_b_max | psi < psi_b_min) {
          
          fprintf(z_map, "%lf, %lf, %lf, %lf, %lf, %lf\n", star.metric.mu[m], phi, nothing, nothing, nothing, nothing);
          fflush(z_map);
          continue;
        }

        int current_p = p;
        double b_final = interp(psi_b, b, MDIV, psi, &current_p);

        if (b_final > b_final_max & phi > 0){
          
          b_final_max = b_final;
          phi_max = phi;
          psi_max = psi;
        }
        if (b_final > b_final_min & phi < 0){

          b_final_min = b_final;
          phi_min = phi;
          psi_min = psi;
        }

        double z_os = redshift_OS(&star, incl_deg, b_final, phi, psi, m);

        double dopp_angle = (1-v_dopp(&star, m)/C*cos_xi(&star, incl_deg, b_final, phi, psi, m));

        double z_con = (1+z_00)*dopp_angle-1;

        double z_qk = (1/zqk_arr[m])*dopp_angle-1;

        double z_qk_bad = (1/zqkb_arr[m])*dopp_angle-1;

        fprintf(z_map, "%lf, %lf, %lf, %lf, %lf, %lf\n", star.metric.mu[m], phi, z_os, z_con, z_qk, z_qk_bad);
        fflush(z_map);

      }

      double z_os_max = redshift_OS(&star, incl_deg, b_final_min, phi_min, psi_min, m);
      double z_os_min = redshift_OS(&star, incl_deg, b_final_max, phi_max, psi_max, m);
      double z_con_min = (1+z_00)*(1-v_dopp(&star, m)/C*cos_xi(&star, incl_deg, b_final_min, phi_min, psi_min, m))-1;
      double z_con_max = (1+z_00)*(1-v_dopp(&star, m)/C*cos_xi(&star, incl_deg, b_final_max, phi_max, psi_max, m))-1;
      
      //fprintf(output, "%d %lf %lf  %lf %lf %lf %lf %lf \n",
	      //m, theta, mu, red_rns, red_sch0, red_sch, red_ker, red_quk);
      
      fprintf(z_out, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n",
        star.metric.mu[m], z_0, z_max, z_min, z_sch, z_grav, z_os_0, z_os_max, z_os_min, gam, red_ker, red_quk, red_quk_fit, z_con_min, z_con_max);
      fflush(z_out);
      //if (m==MDIV) printf("first iteration finished\n");
      //else printf("iteration finished\n");
    }

    fclose(z_out);

    // map of southern hemisphere
    for (m=2;m<=MDIV;m++){

      double z_0 = redshift_metric(&star, 0, m);

      double b_max_OS = b_extreme_OS(&star, m);

      double b_step = (2*b_max_OS)/(MDIV-2);
      double phi_step = (2*PI)/(MDIV-2);
      double psi_b_max = psi_integrated(&star, b_max_OS, m);
      double psi_b_min = psi_integrated(&star, -b_max_OS, m);
      for (p=1;p<MDIV;p++){
        
        double b_trial = -b_max_OS+(p-1)*b_step;
        double psi_trial = psi_integrated(&star, b_trial, m);

        b[p] = b_trial;
        psi_b[p] = psi_trial;
      }

      for (p=1;p<MDIV;p++){

        double phi = -PI+(p-1)*phi_step;

        double psi = acos(-cos(incl)*star.metric.mu[m]+sin(incl)*sqrt(1-SQ(star.metric.mu[m]))*cos(phi));

        if (psi > psi_b_max | psi < psi_b_min) {
          
          fprintf(z_map, "%lf, %lf, %lf, %lf, %lf, %lf\n", -star.metric.mu[m], phi, nothing, nothing, nothing, nothing);
          fflush(z_map);
          continue;
        }

        int current_p = p;
        double b_final = interp(psi_b, b, MDIV, psi, &current_p);

        double z_os = redshift_OS(&star, incl_deg, b_final, phi, psi, m);

        double dopp_angle = (1-v_dopp(&star, m)/C*cos_xi(&star, incl_deg, b_final, phi, psi, m));

        double z_con = (1+z_00)*dopp_angle-1;

        double z_qk = (1/zqk_arr[m])*dopp_angle-1;

        double z_qk_bad = (1/zqkb_arr[m])*dopp_angle-1;

        fprintf(z_map, "%lf, %lf, %lf, %lf, %lf, %lf\n", -star.metric.mu[m], phi, z_os, z_con, z_qk, z_qk_bad);
        fflush(z_map);
      }
    }

    fclose(z_map);
    
    free(b);
    free(psi_b);
    free(zqk_arr);
    free(zqkb_arr);
    }

  return 0;
}









