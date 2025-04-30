#include <stdio.h>
#include <string.h> 
#include <math.h>
#include "nrutil.h"
#include "consts.h"
#include "equil_util.h"
#include "struct.h"

/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/* Adapted from Numerical Recipes.                                         */
/***************************************************************************/
void hunt(double xx[], int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	if ( x >= xx[*jlo] &&  x <= xx[*jlo+1]) {
	  // printf("jlo = %d Perfect! \n",*jlo);
	  return;
	}
	else{

	ascnd=(xx[n] > xx[1]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=0;
		jhi=n+1;
	} else {
		inc=1;
		if ((x >= xx[*jlo]) == ascnd) {
			if (*jlo == n) return;
			jhi=(*jlo)+1;
			while ((x >= xx[jhi]) == ascnd) {
				*jlo=jhi;
				inc += inc;
				jhi=(*jlo)+inc;
				if (jhi > n) {
					jhi=n+1;
					break;
				}
			}
		} else {
			if (*jlo == 1) {
				*jlo=0;
				return;
			}
			jhi=(*jlo);
			*jlo -= 1;
			while ((x < xx[*jlo]) == ascnd) {
				jhi=(*jlo);
				inc += inc;
				*jlo=jhi-inc;
				if (*jlo < 1) {
					*jlo=0;
					break;
				}
			}
		}
	}
	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if ((x > xx[jm]) == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
	}
}

/*C*/
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt)
{ 
 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */ 
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,n_nearest_pt);

 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,1),np+1-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 y= (xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        ((xp[k]-xp[k+1])*(xp[k]-xp[k+2])*(xp[k]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/
       ((xp[k+1]-xp[k])*(xp[k+1]-xp[k+2])*(xp[k+1]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/
       ((xp[k+2]-xp[k])*(xp[k+2]-xp[k+1])*(xp[k+2]-xp[k+3]))
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       ((xp[k+3]-xp[k])*(xp[k+3]-xp[k+1])*(xp[k+3]-xp[k+2]));

 return (y);
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_s(double **f,int s, int m)
{ 
 double d_temp;

 switch(s) { 
            case 1    : d_temp=(f[s+1][m]-f[s][m])/DS;
                        break;

            case SDIV   : d_temp=(f[s][m]-f[s-1][m])/DS;
                          break;
      
            default     : d_temp=(f[s+1][m]-f[s-1][m])/(2.0*DS);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_ss(double **f,int s, int m)
{ 
 double d_temp;

 switch(s) { 
/*
            case 1    : d_temp=(-f[s+3][m]+4.0*f[s+2][m]-5.0*f[s+1][m]
                                                     +2.0*f[s][m] )/SQ(DS);
                        break;

            case SDIV   : d_temp=(2.0*f[s][m]-5.0*f[s-1][m]+4.0*f[s-2][m]
                                                            -f[s-3][m] )/SQ(DS);
                          break;
      
            default     : d_temp=(f[s+1][m]-2.0*f[s][m]+f[s-1][m])/SQ(DS);
                          break; 
*/
/* 
           case 1    : d_temp=(f[s][m]-2.0*f[s+1][m]+f[s+2][m])/(2.0*SQ(DS));
                       break;

           case 2    : d_temp=(2.0*f[s-1][m]-3.0*f[s][m]+f[s+2][m])/(4.0*SQ(DS));
                       break;

           case SDIV-1 : d_temp=(2.0*f[s+1][m]-3.0*f[s][m]+f[s-2][m])
                                                                 /(4.0*SQ(DS));
                         break;

           case SDIV   : d_temp=(f[s][m]-2.0*f[s-1][m]+f[s-2][m])/(2.0*SQ(DS));
                          break;
      
           default     : d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break; 
*/
 
           case 1    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case 2    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case 3    : s=4;
                       d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                       break;

           case SDIV-1 : s=SDIV-2;
                         d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break;

           case SDIV   :  s=SDIV-2;
                          d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                          break;
      
           default     : d_temp=(f[s+2][m]-2.0*f[s][m]+f[s-2][m])/(4.0*SQ(DS));
                         break; 


 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. mu                                */ 
/*******************************************************************/
double deriv_m(double **f,int s, int m)
{
 double d_temp;

 switch(m) { 
            case 1    : d_temp=(f[s][m+1]-f[s][m])/DM;
                        break; 

            case MDIV   : d_temp=(f[s][m]-f[s][m-1])/DM;
                          break;
      
            default     : d_temp=(f[s][m+1]-f[s][m-1])/(2.0*DM);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s                                 */ 
/*******************************************************************/
double deriv_mm(double **f,int s, int m)
{ 
 double d_temp;

 switch(m) { 
            case 1    : m=2;
                        d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                        break;

            case MDIV   : m=MDIV-1;
                          d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                          break;

            default     : d_temp=(f[s][m+1]-2.0*f[s][m]+f[s][m-1])/SQ(DM);
                          break; 
 } 
 return d_temp;
}

/*******************************************************************/
/* Returns the derivative w.r.t. s and mu                          */ 
/*******************************************************************/
double deriv_sm(double **f,int s, int m)
{
 double d_temp;

 switch(s) {
     case 1 : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s][m+1]-f[s+1][m]+f[s][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s][m]-f[s+1][m-1]+f[s][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s+1][m+1]-f[s+1][m-1]-f[s][m+1]+f[s][m-1])/
                                                                (2.0*DM*DS);
                }
              }
              break;

     case SDIV : if(m==1) {   
               d_temp=(f[s][m+1]-f[s][m]-f[s-1][m+1]+f[s-1][m])/(DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s][m]-f[s-1][m]-f[s][m-1]+f[s-1][m-1])/(DM*DS);
                }else{         
                   d_temp=(f[s][m+1]-f[s][m-1]-f[s-1][m+1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }
             }
             break;
  
     default : if(m==1) {   
               d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m]+f[s-1][m])/(2.0*DM*DS);
              }else{
                if(m==MDIV) {
                 d_temp=(f[s+1][m]-f[s-1][m]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (2.0*DM*DS);
                }else{         
                  d_temp=(f[s+1][m+1]-f[s-1][m+1]-f[s+1][m-1]+f[s-1][m-1])/
                                                                (4.0*DM*DS);
                }
             }
             break;
     }

  return d_temp;

}


/*******************************************************************/
/* Returns the Legendre polynomial of degree n, evaluated at x.    */
/*******************************************************************/
double legendre( int n, double x )                      /* checked */
{
  int i;           /* counter */

  double p,        /* Legendre polynomial of order n */
         p_1,      /*    "         "      "    "   n-1*/
         p_2;      /*    "         "      "    "   n-2 */

  p_2=1.0;
  p_1=x;

 if(n>=2) { 
  for(i=2;i<=n;i++){
     p=(x*(2.0*i-1.0)*p_1 - (i-1.0)*p_2)/i;
     p_2=p_1;
     p_1=p;
  }
  return p;
 } else { 
    if (n==1) return p_1;
      else return p_2;
   }
}

/*******************************************************************/
/* Returns the associated Legendre polynomial P_l^m(x).            */
/* Adapted from numerical recipes.                                 */
/*******************************************************************/
double plgndr(int l, int m, double x)
{
	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;

	if (m < 0 || m > l || fabs(x) > 1.0)
		printf("Bad arguments in routine PLGNDR");
	pmm=1.0;
	if (m > 0) {
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++) {
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m)
		return pmm;
	else {
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1))
			return pmmp1;
		else {
			for (ll=(m+2);ll<=l;ll++) {
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return pll;
		}
	}
}

/*C*/
/*******************************************************************/
double rtsec_G(double (*func)(double, double), 
               double Gamma_P, 
               double x1, 
               double x2, 
               double xacc, 
               double ee)
{
 int j;
 double fl,f,dx,swap, xl,rts;
 
 fl=(*func)(x1,Gamma_P)-ee;
 f=(*func)(x2,Gamma_P)-ee;

 if(fabs(fl)<fabs(f)) {
   rts=x1;
   xl=x2;
   swap=fl;
   fl=f;
   f=swap;
 } else {
         xl=x1;
         rts=x2;
        }

 
 for(j=1;j<=MAXIT;j++) {
    dx=(xl-rts)*f/(f-fl);
    xl=rts;
    fl=f;
    rts += dx;
    f=(*func)(rts,Gamma_P)-ee;

    if(fabs(dx)<xacc||f==0.0) return rts;
  }
 
 printf("Maximum number of iterations exceeded in rtsec");  
 return 0.0;
}


/**************************************************************************
 * extrapolate is used to find the surface of the star; s_surf.
 *************************************************************************/

double extrapolate(double *s_edge, double *h_edge, double h_end)
{
  double s_end;
  //if(h_edge[3]==h_end)
  //  h_end -= DBL_EPSILON;
  s_end = ((h_end-h_edge[1])*(h_end-h_edge[2])*(h_end-h_edge[3])*s_edge[0])/
    ((h_edge[0]-h_edge[1])*(h_edge[0]-h_edge[2])*(h_edge[0]-h_edge[3]))
        + ((h_end-h_edge[0])*(h_end-h_edge[2])*(h_end-h_edge[3])*s_edge[1])/
    ((h_edge[1]-h_edge[0])*(h_edge[1]-h_edge[2])*(h_edge[1]-h_edge[3]))
        + ((h_end-h_edge[0])*(h_end-h_edge[1])*(h_end-h_edge[3])*s_edge[2])/
    ((h_edge[2]-h_edge[0])*(h_edge[2]-h_edge[1])*(h_edge[2]-h_edge[3]))
        + ((h_end-h_edge[0])*(h_end-h_edge[1])*(h_end-h_edge[2])*s_edge[3])/
    ((h_edge[3]-h_edge[0])*(h_edge[3]-h_edge[1])*(h_edge[3]-h_edge[2]));
  //printf("s_end = %6.5e = %6.5e  s_edge[3]=%6.5e",
  // ((h_end-h_edge[0])*(h_end-h_edge[1])*(h_end-h_edge[2])*s_edge[3])/
  // ((h_edge[3]-h_edge[0])*(h_edge[3]-h_edge[1])*(h_edge[3]-h_edge[2])),
  //  s_end, s_edge[3]
  //  );

  return s_end;
}

/**************************************************************************
 * interpolate is used to find gama, rho at the surface of the star.
 *************************************************************************/

double interpolate(double *xp, double *yp, double xb)
{
  double yb;
  yb= (xb-xp[1])*(xb-xp[2])*(xb-xp[3])*yp[0]/
        ((xp[0]-xp[1])*(xp[0]-xp[2])*(xp[0]-xp[3]))

    + (xb-xp[0])*(xb-xp[2])*(xb-xp[3])*yp[1]/
       ((xp[1]-xp[0])*(xp[1]-xp[2])*(xp[1]-xp[3]))

    + (xb-xp[0])*(xb-xp[1])*(xb-xp[3])*yp[2]/
       ((xp[2]-xp[0])*(xp[2]-xp[1])*(xp[2]-xp[3]))

    + (xb-xp[0])*(xb-xp[1])*(xb-xp[2])*yp[3]/
       ((xp[3]-xp[0])*(xp[3]-xp[1])*(xp[3]-xp[2]));

  return (yb);
}

// Functions to compute different quantities related to redshift

double r_surf_sch(NeutronStar *star, int mu_i){
   // Schwarzschild radius at a certain mu gridpoint

   double r_sch = star->r_surf[mu_i]*star->Mass*G/SQ(C);

   return r_sch;
}

double r_surf_iso(NeutronStar *star, int mu_i){
   // Isometric radius

   double r_sch = star->r_is_surf[mu_i]*star->Mass*G/SQ(C);

   return r_sch;
}

double b_extreme_metric(NeutronStar *star, int maxmin, int mu_i){
   // b_min (-1) or b_max (+1)

  double term = exp(-star->metric_surf.rho_surf[mu_i])*r_surf_iso(star, mu_i)*sqrt(1-SQ(star->metric.mu[mu_i]));
  double b = maxmin*term/(1+maxmin*star->metric_surf.omega_surf[mu_i]*term/sqrt(KAPPA));

  return b;
}

double v_z(NeutronStar *star, int mu_i){
   // ZAM speed

   double v_z = (star->Omega-star->metric_surf.omega_surf[mu_i]*C/sqrt(KAPPA))*exp(-star->metric_surf.rho_surf[mu_i])*r_surf_iso(star, mu_i)*sqrt(1-SQ(star->metric.mu[mu_i]));

   return v_z;
}

double v_dopp(NeutronStar *star, int mu_i){
   // Doppler speed

   double v = star->Omega*r_surf_sch(star, mu_i)/sqrt(1-(2*star->Mass*G)/(SQ(C)*r_surf_sch(star, mu_i)))*sqrt(1-SQ(star->metric.mu[mu_i]));

   return v;
}

double gamma(NeutronStar *star, int mu_i){
   // Lorentz gamma factor

   double gam = 1/sqrt(1-SQ(v_dopp(star, mu_i)/C));

   return gam;
}

double redshift_metric(NeutronStar *star, double b_z, int mu_i){
   // RNS redshift with b_z

  double z = exp(-(star->metric_surf.gama_surf[mu_i]+star->metric_surf.rho_surf[mu_i])/2)*(1-star->Omega*b_z/C)/(sqrt(1-SQ(v_z(star, mu_i)/C)))-1;
  
  return z;
}

double redshift_OS(NeutronStar *star, double incl_deg, double b, double phi, double psi, int mu_i){
   // Oblate Schwarzschild redshift
   
   double z_os = 1/(sqrt(1-(2*star->Mass*G)/(SQ(C)*r_surf_sch(star, mu_i))))*gamma(star, mu_i)*(1-v_dopp(star, mu_i)/C*cos_xi(star, incl_deg, b, phi, psi, mu_i))-1;
   
   return z_os;
}

double redshift_schwarzschild(NeutronStar *star, int mu_i){
   // Schwarzschild redshift

   double z_sch = 1/(sqrt(1-(2*star->Mass*G)/(SQ(C)*r_surf_sch(star, mu_i))))-1;
   
   return z_sch;
}

double redshift_grav_metric(NeutronStar *star, int mu_i){

   double z_grav = exp(-(star->metric_surf.gama_surf[mu_i]+star->metric_surf.rho_surf[mu_i])/2)-1;

   return z_grav;
}

double redshift_rotation_metric(NeutronStar *star, double b, int mu_i){

   double z_dopp = (1-star->Omega*b/C)/(sqrt(1-SQ(v_z(star, mu_i)/C)))-1;

   return z_dopp;
}

// Angles and parameters for light curves in Oblate Schwarzshild

double b_extreme_OS(NeutronStar *star, int mu_i){

   return r_surf_sch(star, mu_i)/sqrt(1-2*star->Mass*G/(SQ(C)*r_surf_sch(star, mu_i)));
}

double sin_alpha(NeutronStar *star, double b, int mu_i){

   double sin_alpha = b/r_surf_sch(star, mu_i)*sqrt(1-(2*star->Mass*G)/(SQ(C)*r_surf_sch(star, mu_i)));

   return sin_alpha;
}

double psi_integrated(NeutronStar *star, double b, int mu_i){

   double R = r_surf_sch(star, 1);
   double b_hat = b/R;
   double M_R = star->Mass*G/(SQ(C)*R);

   int N = 10000;
   double du = 1/(double)N;
   double u_array[10000];
   for (int i=0; i<N; i++){
      u_array[i] = (i+1/2)*du;
   }

   double integral = (double)0.5 * du/sqrt(1-SQ(b_hat*u_array[0])*(1-2*M_R*u_array[0]));
   for (int i=1; i<(N-1); i++){
      integral += du/sqrt(1-SQ(b_hat*u_array[i])*(1-2*M_R*u_array[i]));
   }
   integral += (double)0.5 * du/sqrt(1-SQ(b_hat*u_array[N-1])*(1-2*M_R*u_array[N-1]));

   return b_hat*integral;
}

double cos_xi(NeutronStar *star, double incl_deg, double b, double phi, double psi, int mu_i){

   if (psi==0) return 0;

   double incl = incl_deg*PI/180;

   return -sin_alpha(star, b, mu_i)*sin(incl)*sin(phi)/sin(psi);
}


double quad_fit(NeutronStar *star){
   // Quadrupole moment 3x3 fit

   double M_R = star->Mass*G/(SQ(C)*r_surf_sch(star, 1));
   double spn = star->Omega*sqrt(r_surf_sch(star, 1)*SQ(r_surf_sch(star, 1))/(star->Mass*G));

   double coefs[3][3] = {
      {  1.27732812e-02, -6.73001462e-03,  6.10940044e-04 },
      {  8.89623481e-01, -1.17221274e-01, -1.20949781e-01 },
      {  1.42132492e-01, -6.09362383e-01,  1.81638146e-01 }
  };

   double quad = 0;

   for (int i=0; i<3; i++){

      for (int j=0; j<3; j++){

         quad += coefs[i][j] * pow(M_R, (double)(-j)) * pow(spn, (double)(2*i));
      }
   }

   return quad;
}
/*****************/
