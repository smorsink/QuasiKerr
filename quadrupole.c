/************************************************************************** 
*                            quadrupole.c                                        
*                                                                         
*       The routines in here are used to find and fit the quardrupole
*       moment for the stellar models. Also includes routines to set up
*       the Boyer-Lindquist metric.                                
*                                                                         
**************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "nrutil.h"
#include "consts.h"
#include "struct.h"
#include "interp.h"
#include "quadrupole.h"
#include "surface.h"

double quadrupole(EOS *eos, NeutronStar *star){
 
  double j;

  // Dimensionless angular momentum
  j = star->ang_mom/pow(star->Mass,2)*C/G;
  star->j = j;

  // Dimensionless Quadrupole Moment (Laarakkers and Poisson)
  star->q = (pow(C, 4.0)*star->Quad)/(SQ(G)*pow(star->Mass, 3.0));

  star->a_factor = star->q/(j*j);

  // Epsilon Factor used in QuasiKerr Metric
  star->epsilon = -star->q - j*j;


  // Quantities used in computation of Coordinate invariant Quadrupole moment.
  star->Btilde *=   pow( star->r_e*C*C/(G*star->Mass),2) * KAPPA;
  star->zeta = star->Btilde + 0.25;
  star->m_scale  = star->q - (4.0/3.0)*star->zeta;
  star->m_2 = star->m_scale * (SQ(G)*pow(star->Mass, 3.0))/pow(C, 4.0) ;

  // star->epsilon = -star->m_scale -j*j;
  
  return star->q;
}

int BLSetUp(NeutronStar *star){
  // Defines the Boyer-Lindquist coordinates used in the Kerr metric.
  // star contains the values describing the RNS spacetime
 
  int i,j,k;

  double *r_BL;
  double *r_iso;
  double *r_S;
  double *driso_dBL;
  double *ds_diso;
  double *g_BL;
  double ***Gamma_BL;

  double mass, a_kerr, r_e, s_gp;

  /* Allocate Memory */

  r_BL = dvector(1,SDIV); // Value of BL radial coordinate at each s grid point in units scaled by the star's mass
  r_iso = dvector(1,SDIV);
  r_S = dvector(1,SDIV);
  driso_dBL = dvector(1,SDIV);
  ds_diso = dvector(1,SDIV);
  g_BL = dvector(0,4);    // Values of BL metric at one value of r_BL
  Gamma_BL = d3tensor(0,3,0,3,0,3);  // Christoffel symbols at one value of r_BL
  
   // mass is the dimensionless mass
   mass = star->Mass * G/(C*C*sqrt(KAPPA));
   a_kerr = star->j; // a = j*mass, but mass=1 in scaled units.
   r_e = star->r_e;

   //  printf("BLSETUP: dimensionless mass = %g  a=%g \n", mass, a_kerr);

   //initialize christoffel symbols
   for (i=0;i<=3;i++)
     for(j=0;j<=3;j++)
       for(k=0;k<=3;k++)
	 Gamma_BL[i][j][k] = 0.0;


   // compute the values of r_BL
   for (i=1; i<=SDIV; i++){

     s_gp = star->metric.s_gp[i];
     r_iso[i] = r_e * s_gp/(1.0-s_gp);
     // convert to "mass" units
     r_iso[i] *= 1.0/mass;

     // Schwarzschild radial coordinate
     r_S[i] = r_iso[i] * exp( 0.5*(star->metric.gama[i][1]-star->metric.rho[i][1]));

     ds_diso[i] = mass/r_e * pow((1.0-s_gp),2);


     if (r_iso[i]==0.0)
       r_BL[i] = 0.0;
     else
       r_BL[i] = r_iso[i] * ( 1.0 + 0.5*(1.0+a_kerr)/r_iso[i]) * ( 1.0 + 0.5*(1.0-a_kerr)/r_iso[i]);

     driso_dBL[i] = 4.0*SQ(r_iso[i])/( 4.0*SQ(r_iso[i]) - (1.0-a_kerr)*(1.0+a_kerr));
   }

  /* Store values in the NeutronStar structure */

  star->blmetric.r_BL = r_BL;
  star->blmetric.r_iso = r_iso;
  star->blmetric.r_S = r_S;
  star->blmetric.driso_dBL = driso_dBL;
  star->blmetric.ds_diso = ds_diso;
  star->blmetric.g_BL = g_BL;
  star->blmetric.Gamma_BL = Gamma_BL;

  return 0;

}



int BLComputeMetric(NeutronStar *star, double rbl, double theta){
  // given a value of Boyer-Lindquist radial and angular coordinates, compute the metric and Christoffel symbols
  // Transforms the RNS metric into the BL coordinate system

  double mu, sintheta;
  double ss;

  double alpha, alpha_s, alpha_m;
  double rho, rho_s, rho_m;
  double gama, gama_s, gama_m;
  double omega, omega_s, omega_m;

  double egama, erho, ealpha, riso, rs, driso_dBL, ds_diso;

  double gtt, gtp, gpp, grr, gthth, det, X;
  double invgrr, invgthth;
  double gprp, gprt, gtrt, gpthp, gptht,gttht, gthrth, grrr, gththth,grthr;
  double a_kerr, mass;

  int snearestpt;
  int first;

  a_kerr = star->j;
  //printf("BLComputeMetric: a_kerr = %g \n", a_kerr);
  mass = star->Mass * G/(C*C*sqrt(KAPPA)); // Dimensionless mass

  if (fabs(theta-0.5*PI) < DBL_EPSILON) {
    mu = 0.0;
    sintheta = 1.0;
  }
  else{
    mu = cos(theta);
    sintheta = sqrt((1.0-mu)*(1.0+mu));
  }

  // rbl = Boyer-Lindquist radial coordinate
  // Given rbl, compute value of s
  // interpolation routines can be found in the file interp.c

  first = (SDIV-1)/2 -1;
  snearestpt=first;

  ss = interp1(star->blmetric.r_BL, star->metric.s_gp, first, SDIV, rbl, &snearestpt);
  riso = interp1(star->blmetric.r_BL, star->blmetric.r_iso, first, SDIV, rbl, &snearestpt);
  rs = interp1(star->blmetric.r_BL, star->blmetric.r_S, first, SDIV, rbl, &snearestpt);
  driso_dBL = interp1(star->blmetric.r_BL, star->blmetric.driso_dBL, first, SDIV, rbl, &snearestpt);
  ds_diso = interp1(star->blmetric.r_BL, star->blmetric.ds_diso, first, SDIV, rbl, &snearestpt);
  
  // Two-dimensional interpolation of the CST metric potentials using Bicubic interpolation. See Numerical Recipes.
  // alpha = value of alpha at the desired Boyer-Linquist coordinate r and theta
  // alpha_s = partial derivative of alpha wrt "s" coordinate. See Cook, Shapiro and Teukolsky
  // alpha_m = partial derivative of alpha wrt "mu=cos(theta)".
  
  pot_interp2( star->metric.alpha, star->metric.alpha_s, star->metric.alpha_mu, star->metric.alpha_ms,
	       ss, mu, &alpha, &alpha_s, &alpha_m);

  pot_interp2( star->metric.rho, star->metric.rho_s, star->metric.rho_mu, star->metric.rho_ms,
	       ss, mu, &rho, &rho_s, &rho_m);

  pot_interp2( star->metric.gama, star->metric.gama_s, star->metric.gama_mu, star->metric.gama_ms,
	       ss, mu, &gama, &gama_s, &gama_m);

  pot_interp2( star->metric.omega, star->metric.omega_s, star->metric.omega_mu, star->metric.omega_ms,
	       ss, mu, &omega, &omega_s, &omega_m);

  // scale omega so that it is measured in mass units

  omega *= mass;
  omega_s *= mass;
  omega_m *= mass;


  // Metric Components in CST Metric (r means r_iso)

  egama = exp(gama);
  erho = exp(rho);
  ealpha = exp(alpha);

  X = pow( riso * sintheta,2) * egama/erho;

  // gtt = -exp(gama+rho) + ( omega r sin(theta) )^2 exp(gama-rho)
  gtt = - egama*erho + pow( omega,2) * X;

  // gtphi = - omega * (r sin(theta))^2 exp(gama-rho)

  gtp = - omega * X;

  // gphiphi = ( r sin(theta) )^2 exp(gama-rho)

  gpp = X;

  // grr (iso) = exp(2 alpha)

  grr = pow( ealpha,2);

  // gthth = r^2 exp(2 alpha)

  gthth = pow(riso,2)*grr;

  star->blmetric.g_BL[0] = gtt;
  star->blmetric.g_BL[1] = grr * pow(driso_dBL,2);
  star->blmetric.g_BL[2] = gthth;
  star->blmetric.g_BL[3] = gpp;
  star->blmetric.g_BL[4] = gtp;

  det = - pow( egama * riso * sintheta,2);

  //contravariant metric components (BL coordinates)
  invgrr = 1.0/star->blmetric.g_BL[1];
  invgthth = 1.0/gthth;

  // lower christoffel symbols

  // gprp = 0.5 * d(gpp)/dr = 0.5 * dX/dr
  gprp = 0.5 * driso_dBL * X * ( 2.0/riso + ds_diso * (gama_s - rho_s));

  // gprt = 0.5 * d(gtp)/dr
  gprt = 0.5 * (- omega * 2.0*gprp - X * driso_dBL*ds_diso * omega_s);

  // gtrt = 0.5*d(gtt)/dr
  gtrt = 0.5 * ( driso_dBL * ds_diso * (-(gama_s + rho_s) * egama*erho + 2.0*omega*omega_s*X)
		 + 2.0*pow(omega,2) * gprp);
		 
  // theta derivatives

  //  gpthp = 0.5 * d(gpp)/d(theta) = 0.5 * dX/dtheta
  if ( mu==1.0)
    gpthp = 0.0;
  else
    gpthp = 0.5 * X * ( 2.0*mu/sintheta - sintheta*(gama_m - rho_m));

  // gptht = 0.5 * d(gpt)/d(theta)
  gptht = 0.5 * ( sintheta*omega_m*X - 2.0*omega*gpthp);

  gttht = 0.5 * ( sintheta * ( (gama_m+rho_m)*egama*erho - 2.0*omega*omega_m*X)
		  + 2.0*pow(omega,2) * gpthp);
  
  // radial derivatives

  // gthrth = 0.5 * d(gthth)/dr
  gthrth = gthth * driso_dBL * ( 1.0/riso + ds_diso * alpha_s);

  // grrr = 0.5 * d(grr)/dr
  grrr = pow(ealpha,2) * pow(driso_dBL,3) * ( ds_diso * alpha_s 
					      - 2.0/riso * (1.0-a_kerr)*(1.0+a_kerr)/( 4*pow(riso,2) - (1.0-a_kerr)*(1.0+a_kerr)));

  // theta derivatives
  gththth = - gthth * sintheta * alpha_m;

  grthr = - pow(ealpha * driso_dBL, 2) * sintheta * alpha_m;

  // Computation of upper Christoffel symbols. 0=t, 1=r, 2=theta, 3=phi
  // First index is the "upper" index. Second and Third indices are the lower indices

  star->blmetric.Gamma_BL[0][1][0] = 0.5 * driso_dBL * ds_diso * ( gama_s + rho_s - omega * omega_s * pow(riso*sintheta/erho,2));
  star->blmetric.Gamma_BL[0][0][1] = star->blmetric.Gamma_BL[0][1][0];

  star->blmetric.Gamma_BL[0][2][0] = -0.5*sintheta * ( (gama_m+rho_m) - omega*omega_m*pow(riso*sintheta/erho,2));
  star->blmetric.Gamma_BL[0][0][2] =   star->blmetric.Gamma_BL[0][2][0];

  star->blmetric.Gamma_BL[0][3][1] = 0.5 * driso_dBL * ds_diso * omega_s * pow(riso*sintheta/erho,2); 
  star->blmetric.Gamma_BL[0][1][3] =   star->blmetric.Gamma_BL[0][3][1];

  star->blmetric.Gamma_BL[0][3][2] = -0.5 * sintheta * omega_m * pow(riso*sintheta/erho,2); 
  star->blmetric.Gamma_BL[0][2][3] =   star->blmetric.Gamma_BL[0][3][2];

  star->blmetric.Gamma_BL[1][0][0] = - invgrr * gtrt;

  star->blmetric.Gamma_BL[1][3][0] = - invgrr * gprt;
  star->blmetric.Gamma_BL[1][0][3] =   star->blmetric.Gamma_BL[1][3][0]; 

  //star->blmetric.Gamma_BL[1][1][1] = invgrr * grrr;
  star->blmetric.Gamma_BL[1][1][1] = driso_dBL * ( ds_diso * alpha_s 
					      - 2.0/riso * (1.0-a_kerr)*(1.0+a_kerr)/( 4*pow(riso,2) - (1.0-a_kerr)*(1.0+a_kerr)));

  star->blmetric.Gamma_BL[1][2][1] = invgrr * grthr;
  star->blmetric.Gamma_BL[1][1][2] =   star->blmetric.Gamma_BL[1][2][1];

  star->blmetric.Gamma_BL[1][2][2] = - invgrr * gthrth;

  star->blmetric.Gamma_BL[1][3][3] = - invgrr * gprp;



  star->blmetric.Gamma_BL[2][0][0] = - invgthth * gttht;

  star->blmetric.Gamma_BL[2][3][0] = - invgthth * gptht;
  star->blmetric.Gamma_BL[2][0][3] =   star->blmetric.Gamma_BL[2][3][0];

  star->blmetric.Gamma_BL[2][1][1] = - invgthth * grthr;

  star->blmetric.Gamma_BL[2][2][1] = invgthth * gthrth;
  star->blmetric.Gamma_BL[2][1][2] =   star->blmetric.Gamma_BL[2][2][1];

  star->blmetric.Gamma_BL[2][2][2] = -sintheta * alpha_m;

  star->blmetric.Gamma_BL[2][3][3] = - invgthth * gpthp;



  star->blmetric.Gamma_BL[3][1][0] = - omega/riso 
    - 0.5 * driso_dBL * ds_diso * ( omega_s - 2.0*omega*rho_s + omega_s * pow(omega*riso*sintheta/erho,2));				
  star->blmetric.Gamma_BL[3][0][1] =   star->blmetric.Gamma_BL[3][1][0];

  star->blmetric.Gamma_BL[3][2][0] = - omega*mu/sintheta +
    0.5 * sintheta * (omega_m - 2.0*omega*rho_m + omega_m* pow(omega*riso*sintheta/erho,2));
  star->blmetric.Gamma_BL[3][0][2] =   star->blmetric.Gamma_BL[3][2][0];

  star->blmetric.Gamma_BL[3][3][1] =  
    driso_dBL * (  1.0/riso + 
		   0.5 * ds_diso * ( gama_s - rho_s + omega*omega_s * pow(riso*sintheta/erho,2)));
  star->blmetric.Gamma_BL[3][1][3] =   star->blmetric.Gamma_BL[3][3][1];

  star->blmetric.Gamma_BL[3][3][2] = mu/sintheta 
    - 0.5 * sintheta * (gama_m - rho_m - omega*omega_m * pow(riso*sintheta/erho,2));
  star->blmetric.Gamma_BL[3][2][3] =   star->blmetric.Gamma_BL[3][3][2];

  return 0;

}

