/* 
************************************************************
             File: metric.c
   Functions that return the various metric elements
   and Christoffel symbols for the spacetime considered.
   This version incorporates the quasi-Kerr metric of 
   Glampedakis & Babak (2006, Class Quant Grav 23, 4167)
   *
   Code from Michi Baubock
   *
************************************************************
*/

#include <stdio.h>
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h> 

#include "metric.h"
//Constants definition file added by Sharon
#include "consts.h"



void metric(double r, double theta, double gmunu[5], double a_bh, double epsilon)
/* March 25, 2010 (DP)
   Function that returns the metric components of the
   spacetime. The coordinates are 'r' and 'theta', while 
   'phi' does not enter because of azimuthal symmetric. The
   properties of the black hole are encoded in the global 
   variables a_bh and epsilon. The non-zero metric components 
   are returned in the array gmunu[5] ordered as 
      gmunu[0]=g_tt
      gmunu[1]=g_rr
      gmunu[2]=g_thth
      gmunu[3]=g_phiphi
      gmunu[4]=g_phit

      SM 20241009: Note that r and a are measured in units of M.
      epsilon = eta x a^2 in the Baubock et al paper.

*/
{
  // Shorthands to be used in the metric elements
  double a2,r2,term_2M_r,costheta,sintheta,cos2theta,costheta2;
  double sintheta2,Logterm;
  
  // The following shorthands should be self explanatory
  a2=a_bh*a_bh;
  r2=r*r;
  term_2M_r=-2.0+r;
  costheta=cos(theta);
  sintheta=sin(theta);
  costheta2=costheta*costheta;
  cos2theta=2.0*costheta2-1.0;
  sintheta2=sintheta*sintheta;
  Logterm=log(r/term_2M_r);
  
  // This is g_tt
  gmunu[0]=
    -((a2 + 2.0*r*term_2M_r + a2*cos2theta)/(a2 + 2.0*r2 + a2*cos2theta)) + 
    (5.0*epsilon*(1.0+3.0*cos2theta)*
     (2.0*(2.0+4.0*r-9.0*r2+3.0*r2*r)-3.0*r2*term_2M_r*term_2M_r*Logterm))/
    (32.0*r2);
  
  // This is g_rr
  gmunu[1]=
    (r2 + a2*costheta2)/(a2-2.0*r+r2) - 
    (5.0*epsilon*(1.0-3.0*costheta2)*(2.0*(1.0-r)*(2.0+6.0*r-3.0*r2)-3.0*r2*term_2M_r*term_2M_r*Logterm))/
    (16.0*term_2M_r*term_2M_r);
  
  // This is g_thth
  gmunu[2]=
    r2 + a2*costheta2 - 
    (5.0*r*epsilon*(1.0+3.0*cos2theta)*
     (2.0*(2.0-3.0*r-3.0*r2)+(-6.0*r+3.0*r2*r)*Logterm))/(32.0);
  
  // This is g_phiphi
  gmunu[3]=
    (-5.0*r*epsilon*(1.0+3.0*cos2theta)*
	   (2.0*(2.0-3.0*r-3.0*r2)+(-6.0*r+3.0*r2*r)*Logterm)*sintheta2)
    /(32.0) + 
    (4.0*(r2+a2*costheta2)*sintheta2*
     ((a2+r2)*(a2+r2)-a2*(a2 + r*term_2M_r)*sintheta2))/
    ( (a2+2.0*r2+a2*cos2theta)*(a2+2.0*r2+a2*cos2theta));
  
  // This is g_tphi
  gmunu[4]=
    (-4.0*a_bh*r*sintheta2)/(a2+2.0*r2+a2*cos2theta);
  
  return;
  // all done
}

void gammas(double r, double theta, double Gamma[4][4][4], double a_bh, double epsilon)
/* March 30, 2010 (DP), based on the subroutine by TJ.
   Function that returns the Christoffel symbols for the
   spacetime. The coordinates are 'r' and 'theta', while 
   'phi' does not enter because of azimuthal symmetric. The
   properties of the black hole are encoded in the global
   variable a_bh and epsilon. The Christoffel symbols are 
   returned in the array Gamma[4][4][4]. The ordering of 
   the indices for the symbols is
   0=t, 1=r, 2=theta, 3=phi */
{
  // Shorthands to be used for the Christoffel symbols
  double sintheta, costheta, sintheta2, costheta2,sincos;
  double r2, a2, Delta, Sigma, Sigma2, mSigma;
  double ra2, r3, term_2M_r, term2, Logterm, c2s, poly5, poly6, poly7, poly8;

  // The following shorthands should be self explanatory
  // ... about angles

  // Small section added by Sharon so that values on equatorial plane are correct.
  if (fabs(theta-0.5*PI) < DBL_EPSILON) {
    costheta = 0.0;
    sintheta = 1.0;
  }
  else{
    costheta = cos(theta);
    sintheta = sin(theta);
  }

  sintheta2 = sintheta*sintheta;
  costheta2 = costheta*costheta;
  sincos = sintheta*costheta;
  c2s = 1.0 + 3.0*cos(2.0*theta);
  // ... about the radial coordinates
  r2 = r*r;  
  a2 = a_bh*a_bh;
  ra2 = r2 + a2;
  r3 = r2*r;
  Delta = r2 + a2 - 2.0*r;
  Sigma = r2 + a2 * costheta2;
  term_2M_r=-2.0+r;
  term2=term_2M_r*term_2M_r;
  Logterm=log(r/term_2M_r);
  Sigma2 = Sigma*Sigma;
  mSigma = r2-a2*costheta2;
  // ... and a few polynomials that appear in the expressions
  poly5 = 5.0*epsilon*c2s * 
    (2.0 * (2.0* - 2.0*r + 13.0*r2 - 12.0*r3 + 3.0*r2*r2) + 
     3.0*(1.0 - r)*r2*term2*Logterm) / (32.0*r2*term2);
  poly6 = 5.0*epsilon*c2s * (2.0 * (14.0 - 3.0*r - 3.0*r2) + 
			   3.0*(4.0 - 6.0*r + r3)*Logterm) / (32.0);
  poly7 = 15.0*epsilon * (2.0 * (2.0 - 3.0*r - 3.0*r2) + 
		       3.0*r * (r2 - 2.0)*Logterm)*sincos / (16.0*r);
  poly8 = 5.0*epsilon*c2s * 
    (2.0 * (2.0 + 2.0*r + 3.0*r2 - 3.0*r3) + 3.0*r3*term_2M_r*Logterm) / 
    (32.0*r2*term_2M_r);
  

  Gamma[0][1][0] = ra2*mSigma/ ( Delta*Sigma2 ) - poly5;
  Gamma[0][0][1] = Gamma[0][1][0];    // Symmetry of lower two indices

  Gamma[0][2][0] = -2.0*a2*r*sincos/Sigma2 + 15.0*epsilon * 
    (2.0*(1.0 - r)*(2.0 + 6.0*r - 3.0*r2) - 3.0*r2*term2*Logterm) 
    * sincos / (16.0*r*term_2M_r);
  Gamma[0][0][2] = Gamma[0][2][0];    // Symmetry of lower two indices
  
  


  Gamma[0][3][1] = a_bh*sintheta2 * ( (r2-a2)*mSigma - 4.0*r2*r2 )
    /( Delta*Sigma2 );
  Gamma[0][1][3] = Gamma[0][3][1];    // Symmetry of lower two indices

  Gamma[0][3][2] = 2.0*a_bh*a2*r*sintheta2*sincos /Sigma2;
  Gamma[0][2][3] = Gamma[0][3][2];    // Symmetry of lower two indices

  Gamma[1][0][0] = Delta*mSigma/(Sigma*Sigma2) + 
    5.0*epsilon*c2s * (-2.0 * (6.0 + 6.0*r - 5.0*r2 - 6.0*r3 + 3.0*r2*r2) 
		    + 3.0*r2*term2*(1.0+r)*Logterm) / (32.0*r2*r2);

  Gamma[1][3][0] = -a_bh*sintheta2*mSigma*Delta / (Sigma*Sigma2);
  Gamma[1][0][3] = Gamma[1][3][0];    // Symmetry of lower two indices

  Gamma[1][1][1] = (1.0 - r)/Delta + r/Sigma + poly5;
  
  Gamma[1][2][1] = -a2*sincos/Sigma - 
    15.0*epsilon * (2.0 * (2.0 + 4.0*r - 9.0*r2 + 3.0*r3) - 
		3.0*r2*term2*Logterm)*sincos / (16.0*r*term_2M_r);
  Gamma[1][1][2] = Gamma[1][2][1];    // Symmetry of lower two indices

  Gamma[1][2][2] = -r*Delta/Sigma + poly6;
  
  Gamma[1][3][3] = Delta*sintheta2*
    ( 2.0*r2*ra2 - (r-1.0)*Sigma2 - (2.0*r2+ra2)*Sigma ) 
    / (Sigma*Sigma2) + poly6*sintheta2;

  Gamma[2][0][0] = -2.0*a2*r*sincos/(Sigma*Sigma2) + 
    15.0*epsilon * (2.0 * (2.0 + 4.0*r - 9.0*r2 + 3.0*r3) 
		- 3.0*r2*term2*Logterm)*sincos / (16.0*r2*r2);

  Gamma[2][3][0] = 2.0*a_bh*r*ra2*sincos / (Sigma*Sigma2);
  Gamma[2][0][3] = Gamma[2][3][0];    // Symmetry of lower two indices
  
  Gamma[2][1][1] = a2*sincos/( Delta*Sigma ) + 
    15.0*epsilon*sincos * (2.0*(1. - r) * (2.0 + 6.0*r - 3.0*r2) 
		      - 3.0*r2*term2*Logterm) / (16.0*r2*term2);

  Gamma[2][2][1] = r/Sigma - poly8;
  Gamma[2][1][2] = Gamma[2][2][1];    // Symmetry of lower two indices
  
  Gamma[2][2][2] = -a2*sincos/Sigma + poly7;
  
  Gamma[2][3][3] = - sincos * ( Delta*Sigma2 + 2.0*r*ra2*ra2 ) 
    / (Sigma*Sigma2) - poly7*sintheta2;

  Gamma[3][1][0] = a_bh*mSigma/(Delta*Sigma2);
  Gamma[3][0][1] = Gamma[3][1][0];    // Symmetry of lower two indices
  
  Gamma[3][2][0] = -2.0*a_bh*r*costheta/(sintheta*Sigma2);
  Gamma[3][0][2] = Gamma[3][2][0];    // Symmetry of lower two indices

  Gamma[3][3][1] = ( r*(Sigma-2.0*r)*Sigma - 
		     mSigma*a2*sintheta2 ) / (Delta*Sigma2) - poly8;
  Gamma[3][1][3] = Gamma[3][3][1];    // Symmetry of lower two indices

  Gamma[3][3][2] = (Sigma2+2.0*a2*r*sintheta2) * costheta 
    / (sintheta*Sigma2) + 15.0*epsilon * 
    (2.0* (2.0 - 3.0*r - 3.0*r2) + (3.0*r3 - 6.0*r)*Logterm)*sincos / (16.0*r);
  Gamma[3][2][3] = Gamma[3][3][2];    // Symmetry of lower two indices

  return;
}

