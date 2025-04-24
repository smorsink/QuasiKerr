#include "struct.h"

void hunt(double xx[], int n, double x, int *jlo);

double interp(double xp[], 
              double yp[], 
              int    np ,
              double xb, 
              int    *n_nearest_pt);

double deriv_s(double **f,int s, int m);

double deriv_ss(double **f,int s, int m);

double deriv_m(double **f,int s, int m);

double deriv_mm(double **f,int s, int m);

double deriv_sm(double **f,int s, int m);

double legendre( int n, double x );

double plgndr(int l, int m, double x);
 
double rtsec_G( double (*func)(double, double), 
                double Gamma_P,
                double x1, 
                double x2, 
                double xacc,
                double ee);

double extrapolate(double *s_edge, double *h_edge, double h_end);

double interpolate(double *xp, double *yp, double xb);

double r_surf_sch(NeutronStar *star, int mu_i);

double b_extreme_metric(NeutronStar *star, int maxmin, int mu_i);

double redshift_metric(NeutronStar *star, double b, int mu_i);

double v_z(NeutronStar *star, int mu_i);

double v_dopp(NeutronStar *star, int mu_i);

double gamma(NeutronStar *star, int mu_i);

double redshift_OS(NeutronStar *star, double incl_deg, double b, double phi, double psi, int mu_i);

double redshift_schwarzschild(NeutronStar *star, int mu_i);

double redshift_grav_metric(NeutronStar *star, int mu_i);

double redshift_rotation_metric(NeutronStar *star, double b, int mu_i);

double b_extreme_OS(NeutronStar *star, int mu_i);

double psi_integrated(NeutronStar *star, double b, int mu_i);

double cos_xi(NeutronStar *star, double incl_deg, double b, double phi, double psi, int mu_i);

double quad_fit(NeutronStar *star);