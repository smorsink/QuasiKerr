/***************************************************************************
 * interp.c
 *
 * This file contains a collection of useful interpolation routines.
 *
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "consts.h"
#include "nrutil.h"
#include "equil_util.h" /*get hunt and interp from here*/
#include "interp.h"

/* Do a linear interpolation between two points. */
double interplin(double *xp, double *yp, int np, double xb, int *n_nearest_pt){
  int k;

  hunt(xp,np,xb,n_nearest_pt);
  k=*n_nearest_pt;
  return yp[k] + (xb-xp[k])*(yp[k+1]-yp[k])/(xp[k+1]-xp[k]);
}

/* Polynomial interpolation.  This is based on the interpolation routine in
 * numerical recipes.  If P(x) is the polynomial of order order-1 passing 
 * through the points (xp[i],yp[i]), then polint returns P(xb) and err points
 * to an error estimate.  Note xp and yp are arrays indexed from 1 to n>=order.
 */
double polint(double *xp, double *yp, int order, double xb, double *err){
  int i, m, ns=1;
  double tmp, diff, den, dnum, cnum, yb;
  double *c, *d;

  c = (double *) malloc((order+1)*sizeof(double));
  d = (double *) malloc((order+1)*sizeof(double));
  diff = fabs(xb-xp[1]);
  for (i=1; i<=order; i++){
    if ( (tmp=fabs(xb-xp[i])) < diff ){
      ns=i;
      diff = tmp;
    }
    c[i] = yp[i];
    d[i] = yp[i];
  }

  yb = yp[ns--];
  for (m=1; m<order; m++){
    for (i=1; i<=order-m; i++){
      cnum = xp[i]-xb;
      dnum = xp[i+m]-xb;
      if ( (den=cnum-dnum) == 0.0 ){
        /*Two values of xp are equal:no polynomial passes through the points.*/
        printf("error in polint: xp[%d]==xp[%d]\n", i, i+m);
        exit(1);
      }
      tmp = (c[i+1]-d[i])/den;
      c[i] = cnum*tmp;
      d[i] = dnum*tmp;
    }
    *err = 2*ns<order-m ? c[ns+1] : d[ns--];
    yb += *err;
  }
  free(c);
  free(d);
  return yb;
}

/* Do a 1-dimensional interpolation. */
double interp1(double *xp, double *yp, int first, int np, double xb, int *x_nearest_pt){
  const int order = 4; /*Order of interpolation*/
  int lo;

  double err;
  double ans;

  hunt_surf(xp,first,np,xb,x_nearest_pt);

  lo=IMIN(IMAX(*x_nearest_pt-(order-1)/2, 1), np+1-order);

  ans = polint(&xp[lo-1], &yp[lo-1], order, xb, &err);

  //printf("ans = %g err=%g \n",ans,err);

  return ans;

}


// hunt_surf searches the region outside of the star to find the appropriate value
void hunt_surf(double xx[], int first, int n, double x, int *jlo)
{ 
	int jm,jhi,inc,ascnd;

	if ( x >= xx[*jlo] &&  x <= xx[*jlo+1]) {
	  //printf("jlo = %d Perfect! \n",*jlo);
	  //printf("x = %g  ; xx[jlo] = %g ; xx[jlo+1] = %g \n",
	  //	 x, xx[*jlo], xx[*jlo+1]);
	  return;
	}
	else{

	ascnd=(xx[n] > xx[first]);
	if (*jlo <= 0 || *jlo > n) {
		*jlo=first-1;
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
			if (*jlo == first) {
				*jlo=first-1;
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

	//	printf("x = %g jhi = %d xx[jhi] = %g \n", x, jhi, xx[jhi]);

	while (jhi-(*jlo) != 1) {
		jm=(jhi+(*jlo)) >> 1;
		if ((x > xx[jm]) == ascnd)
			*jlo=jm;
		else
			jhi=jm;
	}
	}
}



/* Do a 2-dimensional interpolation.  */
double interp2(double *xp, double *yp, double **zp, int nx, int ny,
               double xb, double yb, int *x_nearest_pt, int *y_nearest_pt){

  const int max_order = 4;  /*Interpolation in the x-direction is of the order
                              max(max_order, m) where m is the order of 
                              interpolation used in interp.*/
  int i,lo;
  int nnp;
  double y_interp[max_order+1];
  
  hunt(xp, nx, xb, x_nearest_pt);
  lo = IMIN(IMAX(*x_nearest_pt-(max_order-1)/2, 1), nx+1-max_order);
  for (i=1; i<=max_order; i++){
    y_interp[i] = interp(yp, zp[lo+i-1], ny, yb, y_nearest_pt);
  }
  nnp = max_order >> 1;
  return interp(&xp[lo-1], y_interp, max_order, xb, &nnp);
}


/* Do a 1-D interpolation on the first index of a 2-D array. */
double interp_b(double *xp, double **zp, int nx, int y, double xb,
                int *n){
  int i;
  double *z_at_y, val;


  z_at_y = dvector(1, nx);
  for (i=1; i<=nx; i++)
    z_at_y[i] = zp[i][y];
  val =  interp(xp, z_at_y, nx, xb, n);
  free_dvector(z_at_y, 1, nx);
  return val;
}







/* New version of interp2 Do a 2-dimensional interpolation.  */
double interp2A(double *xp, double *yp, double **zp, 
               double xb, double yb, int xlo, int ylo, int max_order){

  int i;
  double y_interp[max_order+1];
  double  ans;
   
  for (i=1; i<=max_order; i++){
    
        y_interp[i] = interpB(yp, zp[xlo+i-1], yb, ylo);
    
  }

 if( xb==xp[xlo] ||  xb==xp[xlo+1] || xb==xp[xlo+2] || xb==xp[xlo+3]) 
    xb += DBL_EPSILON;

  ans = -(xb-xp[xlo+1])*(xb-xp[xlo+2])*(xb-xp[xlo+3])*y_interp[1]/3.0
    +(xb-xp[xlo])*(xb-xp[xlo+2])*(xb-xp[xlo+3])*y_interp[2]
    -(xb-xp[xlo])*(xb-xp[xlo+1])*(xb-xp[xlo+3])*y_interp[3]
    +(xb-xp[xlo])*(xb-xp[xlo+1])*(xb-xp[xlo+2])*y_interp[4]/3.0;
      
  ans *= 0.5*pow(1.0/DS,3);

 return ans;
  
}


/*C*/
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */  
/*************************************************************************/
double interpB(double xp[], 
              double yp[], 
              double    xb ,
              int    k)
{ 
 
 double y;     /* intermediate value */


 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3]) 
    xb += DBL_EPSILON;

 y= -(xb-xp[k+1])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k]/
        (6.0)
 
    +(xb-xp[k])*(xb-xp[k+2])*(xb-xp[k+3])*yp[k+1]/2.0
       
 
    -(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+3])*yp[k+2]/2.0
       
 
    +(xb-xp[k])*(xb-xp[k+1])*(xb-xp[k+2])*yp[k+3]/
       (6.0);

 y *= pow(1.0/DM,3);

 return (y);
}


/* Interpolate to find the value of a potential at r = infinity */

double interp_infty(double **pot, double mu) 
{ 

  int i,k=4,left_mu,lo;
  double ans;     /* intermediate value */

  static double mupts[4],potpts[4];


  if (mu == 0.0) ans= pot[SDIV][1]; /* equator */  
  else
    if (mu == 1.0) ans= pot[SDIV][MDIV]; /* North Pole */
    else{ /* All other points */

      left_mu = mu/DM + 1;

      lo = IMIN(IMAX(left_mu-(k-1)/2, 1), MDIV+1-k);

      for (i=0;i<k;i++){
	mupts[i] = (lo+i-1.0-1.0)/(MDIV-1.0);
	potpts[i] = pot[SDIV][lo+i-1];
      }
      
      k = 0;

      ans = -(mu-mupts[k+1])*(mu-mupts[k+2])*(mu-mupts[k+3])*potpts[k]/
        (6.0)
 
	+(mu-mupts[k])*(mu-mupts[k+2])*(mu-mupts[k+3])*potpts[k+1]/2.0
       
 
	-(mu-mupts[k])*(mu-mupts[k+1])*(mu-mupts[k+3])*potpts[k+2]/2.0
       
 
	+(mu-mupts[k])*(mu-mupts[k+1])*(mu-mupts[k+2])*potpts[k+3]/
       (6.0);
      
      ans *= pow(1.0/DM,3);

    }

  return ans;
 
}





/* Copyright (C) 1987,1988 Numerical Recipes Software -- BCUCOF */

void bcucof(double *y,
	    double *y1,
	    double *y2,
	    double *y12,
	    double d1,
	    double d2,
	    double **c)

{
  static int wt[16][16]=
    {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
     {-3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0},
     {2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0},
     {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
     {0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1},
     {0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1},
     {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
     {9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2},
     {-6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2},
     {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
     {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
     {-6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1},
     {4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1}};

	int l,k,j,i;
	double xx,d1d2,cl[16],x[16];


    
	d1d2=d1*d2;
	for (i=1;i<=4;i++) {
		x[i-1]=y[i];
		x[i+3]=y1[i]*d1;
		x[i+7]=y2[i]*d2;
		x[i+11]=y12[i]*d1d2;
	}
	for (i=0;i<=15;i++) {
		xx=0.0;
		for (k=0;k<=15;k++) xx += wt[i][k]*x[k];
		cl[i]=xx;
	}
	l=0;
	for (i=1;i<=4;i++)
		for (j=1;j<=4;j++) c[i][j]=cl[l++];
}


/* Copyright (C) 1987,1988 Numerical Recipes Software -- BCUINT */

void bcuint(double *y,
	    double *y1,
	    double *y2,
	    double *y12,
	    double x1l,
	    double x1u,
	    double x2l,
	    double x2u,
	    double x1,
	    double x2,
	    double *ansy,
	    double *ansy1,
	    double *ansy2)

{
	int i;
	double t,u,d1,d2,**c;

   
	c=dmatrix(1,4,1,4);
	d1=x1u-x1l;
	d2=x2u-x2l;
	bcucof(y,y1,y2,y12,d1,d2,c);


	if (x1u == x1l || x2u == x2l) nrerror("Bad input in routine BCUINT");
	t=(x1-x1l)/d1;
	u=(x2-x2l)/d2;
	(*ansy)=0.0;
	(*ansy2)=0.0;
	(*ansy1)=0.0;
	for (i=4;i>=1;i--) {
		*ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
		*ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
		*ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];
	}
	*ansy1 /= d1;
	*ansy2 /= d2;
	free_dmatrix(c,1,4,1,4);



}


/* Interpolate a funtion on the equator */

double interp_equat(double **pot, double ss)
{

  int i,k=4,left_s,lo;
  double y;     /* intermediate value */

  static double xp[5],yp[5];

  if (ss == 0.5) left_s = ss/DS;
  else
    left_s = ss/DS + 1;

  lo = IMIN(IMAX(left_s-(k-1)/2, 1), SDIV+1-k);

  for (i=1;i<=k;i++){

    xp[i] = SMAX*(lo+i-1.0-1.0)/(SDIV-1.0);
    yp[i] = pot[lo+i-1][1];

    //printf("i=%d xp[i]=%6.5e yp[i]=%6.5e",i,xp[i],yp[i]);

  }

 if( ss==xp[k] ||  ss==xp[k+1] || ss==xp[k+2] || ss==xp[k+3]) 
    ss += DBL_EPSILON;

 y= -(ss-xp[2])*(ss-xp[3])*(ss-xp[4])*yp[1]/
        (6.0)
 
    +(ss-xp[1])*(ss-xp[3])*(ss-xp[4])*yp[2]/2.0
       
 
    -(ss-xp[1])*(ss-xp[2])*(ss-xp[4])*yp[3]/2.0
       
 
    +(ss-xp[1])*(ss-xp[2])*(ss-xp[3])*yp[4]/
       (6.0);

 y *= pow(1.0/DS,3);

 //printf("y = %6.5e \n",y);

 return (y);


}






/* New version of interp2 Do a 2-dimensional interpolation.  */
void    pot_interp2(		   
		    double **pot, double **pot_s, double **pot_m, double **pot_ms,
		    double ss, double mm, 
		    double *ans, double *ans_s, double *ans_m)
{

  int left_s, left_m;

  double y[5], y1[5], y2[5], y12[5];
  double x1l,x2l,x1u,x2u;

  if (ss == 0.5) left_s = ss/DS;
  else
    left_s = ss/DS + 1;
  left_m = mm/DM + 1;

  if (left_s+1 > SDIV){
    //printf("error in pot_interp!! ss=%g SDIV=%d DS=%g left_s=%d = %d\n",
    //   ss, SDIV, DS, left_s, ss/DS+1);
    left_s = SDIV-1;
  }


  y[1] = pot[left_s][left_m];
  y[2] = pot[left_s+1][left_m];
  y[3] = pot[left_s+1][left_m+1];
  y[4] = pot[left_s][left_m+1];


  y1[1] = pot_s[left_s][left_m];
  y1[2] = pot_s[left_s+1][left_m];
  y1[3] = pot_s[left_s+1][left_m+1];
  y1[4] = pot_s[left_s][left_m+1];

  y2[1] = pot_m[left_s][left_m];
  y2[2] = pot_m[left_s+1][left_m];
  y2[3] = pot_m[left_s+1][left_m+1];
  y2[4] = pot_m[left_s][left_m+1];

  y12[1] = pot_ms[left_s][left_m];
  y12[2] = pot_ms[left_s+1][left_m];
  y12[3] = pot_ms[left_s+1][left_m+1];
  y12[4] = pot_ms[left_s][left_m+1];

  x1l = DS*(left_s-1);
  x1u = x1l+DS;
  
  x2l = DM*(left_m-1);
  x2u = x2l+DM;

  bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,ss,mm,ans,ans_s,ans_m);



}


/* Polynomial Extrapolation.  This is based on the interpolation routine in
 * numerical recipes.  If P(x) is the polynomial of order order-1 passing 
 * through the points (xp[i],yp[i]), then polint returns P(xb) and err points
 * to an error estimate.  Note xp and yp are arrays indexed from 1 to n>=order.
 */
double extrap(double *xp, double *yp, int order, double xb, double *err){
  int i, m, ns;
  double tmp, den, dnum, cnum, yb;
  double *c, *d;

  c = (double *) malloc((order+1)*sizeof(double));
  d = (double *) malloc((order+1)*sizeof(double));

  ns = order;


  for (i=1; i<=order; i++){
    c[i] = yp[i];
    d[i] = yp[i];
  }

  yb = yp[ns--];
  for (m=1; m<order; m++){
    for (i=1; i<=order-m; i++){
      cnum = xp[i]-xb;
      dnum = xp[i+m]-xb;
      if ( (den=cnum-dnum) == 0.0 ){
        /*Two values of xp are equal:no polynomial passes through the points.*/
        printf("error in polint: xp[%d]==xp[%d]\n", i, i+m);
        exit(1);
      }
      tmp = (c[i+1]-d[i])/den;
      c[i] = cnum*tmp;
      d[i] = dnum*tmp;
    }
    *err = 2*ns<order-m ? c[ns+1] : d[ns--];
    yb += *err;
  }
  free(c);
  free(d);
  return yb;
}
