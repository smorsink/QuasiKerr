#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <boost/math/tools/roots.hpp>
#include <spline.h/spline.h>

using namespace std;





// constants for user input

const string B_EOS = "NL3wrL55_500"; //"APR" //"SLy5"; //"NL3wrL55"; // baryonic eos

// MeVfm constants
//const double m_phi_input_MeVfm = 1.0e-3; // DM self-interaction mediator mass
const double m_D_input_MeVfm = 1.0e3;            // DM particle mass
const double epsilon_B_c_input_MeVfm = 365.295;  // baryonic central energy density
const double epsilon_D_c_input_MeVfm = 0.0641385; // DM central energy density

// dimensionless constants
const int B_EOS_CompOSE_flag = 1;        // 1 if CompOSE, 0 if rns
//const double f_D_target_input = 0.05; // desired DM fraction of DANS
//const double g_input = 1.0e-4; // DM and mediator coupling strength
const double y_input = 1000.0; //0.0; //1.0e3;       // DM self-interaction strength


// constants not to be changed
//***********CHECK CONVERGENCE BY CHANGING FOLLOWING VALUES*********
// - r_TOV_accurateregion_max_km
// - dr_TOV_min_km
// - dr_TOV_max_km
// - r_powerseries_max_factor
// - dr_TOV_accurateregion_max_factor
// - surf_dr_TOV_shrink_factor
// - x_max
// - mass_rootfind_accuracy
// - P_B_c_rootfind_accuracy
// - x_rootfind_accuracy
// - RKF45_tolerance_factor

// km constants
const double r_TOV_accurateregion_max_km = 25.0; // upper limit for high resolution region for solving TOV equations
const double dr_TOV_min_km = 1.0e-3;             // minimum grid spacing for solving TOV equations
const double dr_TOV_max_km = 1.0e2;             // minimum grid spacing for solving TOV equations

// MeVfm constants
const double GeV_MeVfm = 1.0e3;                           // GeV
//const double m_phi_MeVfm = m_phi_input_MeVfm; // DM self-interaction mediator mass
const double m_D_MeVfm = m_D_input_MeVfm;                 // DM particle mass
const double epsilon_B_c_MeVfm = epsilon_B_c_input_MeVfm; // baryonic central energy density
const double epsilon_D_c_MeVfm = epsilon_D_c_input_MeVfm; // DM central energy density

// cgs constants
const double MeV_cgs = 1.60218e-6;       // MeV
const double fm_cgs = 1.0e-13;           // fm
const double km_cgs = 1.0e5;             // km
const double solarmass_cgs = 1.989e33;   // solar mass
const double c_cgs = 2.99792458e10;      // speed of light
const double G_cgs = 6.6743e-8;          // Newton's gravitational constant
const double hbar_cgs = 1.054571817e-27; // reduced Planck's constant
const double h_B_surf_cgs = 1.0;         // baryonic enthalpy at surface

// dimensionless constants
const double r_powerseries_max_factor = 4.0e0;                   // factor of dr_TOV_min to set maximum radius up to which to use powerseries
const double dr_TOV_accurateregion_max_factor = 5.0e2;           // factor of dr_TOV_min to set maximum grid spacing in high resolution region
const double dr_TOV_max_factor = 1.0e-2;                         // factor of the order of magnitude of radius to have maximum grid spacing
const double dr_TOV_surf_max_factor = 5.0e-4;                    // factor of dr_TOV_min to set maximum grid spacing in high resolution region
const double x_min = 0.0;                                        // minimum possible value of x
const double x_max = 1.0e6;                                      // maximum value of x
//const double g = g_input;                           // DM and mediator coupling strength
//const double y = g*m_D_MeVfm/(m_phi_MeVfm*sqrt(2)); // DM self-interaction strength
const double y = y_input;                                        // DM self-interaction strength;
const double r2_c = 0.0;                                          // central value of r
const double epsilon_B_c_zero = 0.0;                             // zero baryonic central energy density
const double M_B_c = 0.0;                                        // central baryonic mass
const double N_B_c = 0.0;                                        // central number of baryons
const double epsilon_B_surf_zero = 0.0;                               // baryonic energy density at surface
const double P_B_surf_zero = 0.0;                                // zero baryonic pressure at surface
const double h_B_surf_zero = 0.0;                                     // baryon number density at surface
const double n_B_surf_zero = 0.0;                                     // baryon number density at surface
const double epsilon_D_c_zero = 0.0;                             // zero DM central energy density
const double M_D_c = 0.0;                                        // central DM mass
const double N_D_c = 0.0;                                        // central number of DM particles
const double epsilon_D_surf_zero = 0.0;                               // DM energy density at surface
const double P_D_surf_zero = 0.0;                                // zero DM pressure at surface
const double h_D_surf_zero = 0.0;                                     // DM particle number density at surface
const double n_D_surf_zero = 0.0;                                     // DM particle number density at surface
const double f_c = 1.0;                                          // central value of rr metric
const double R_D_zero = 0.0;                                     // minimum value of R_D
const double M_halo_over_R_D_zero = 0.0;                         // minimum value of M_halo/R_D
const double M_halo_over_R_D_minus_R_B_zero = 0.0;
const double M_halo_over_M_D_of_R_D_zero = 0.0;                  // minimum value of M_halo/R_D
const double s_denom_zero = 0.0;                                 // zero s denominator
const double P_B_c_rootfind_accuracy = 1.0e-10;                  // baryonic central pressure root finding accuracy
const double x_rootfind_accuracy = 1.0e-15;                      // x root finding accuracy
const double RKF45_tolerance_factor = 1.0e-10;                    // RKF45 tolerance factor
const double RKF45_tolerance_factor_surf = 1.0e-3;               // RKF45 tolerance factor near surface
const int R_D_index_zero = 0;                                    // surface index for R_D = 0


// unit conversion of constants to cgs
//const double m_phi_cgs = m_phi_MeVfm*MeV_cgs/pow(c_cgs, 2.0); // DM self-interaction mediator mass
const double m_D_cgs = m_D_MeVfm*MeV_cgs/pow(c_cgs, 2.0);                                    // DM particle mass
const double r_TOV_accurateregion_max_cgs = r_TOV_accurateregion_max_km*km_cgs;              // upper limit for high resolution region for solving TOV equations
const double dr_TOV_min_cgs = dr_TOV_min_km*km_cgs;                                          // minimum grid spacing for solving TOV equations
const double dr_TOV_max_cgs = dr_TOV_max_km*km_cgs;                                          // minimum grid spacing for solving TOV equations
const double epsilon_B_c_cgs = epsilon_B_c_MeVfm*MeV_cgs/(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)); //epsilon_B_c_MeVfm*MeV_cgs/(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)); // baryoninc central energy density
const double epsilon_D_c_cgs = epsilon_D_c_MeVfm*MeV_cgs/(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)); // epsilon_D_c_MeVfm*MeV_cgs/(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)); // DM central energy density


// declare constants for unit conversion to dimensionless
double hbar;                      // reduced Planck's constant
//const double m_phi; // DM self-interaction mediator mass
double r_TOV_accurateregion_max;  // upper limit for high resolution region for solving TOV equations
double dr_TOV_min;                // minimum grid spacing for solving TOV equations
double dr_TOV_max;                // minimum grid spacing for solving TOV equations
double r_powerseries_max;         // end point of power series
double dr_TOV_accurateregion_max; // maximum grid spacing in high resolution region for solving TOV equations
double dh;
double P_B_surf;                  // baryonic pressure at surface
double h_B_surf;                  // baryonic enthalpy at surface

// declare RKF45 tolerances
double RKF45_tolerance_B; // baryonic RKF45 tolerance
double RKF45_tolerance_D; // DM RKF45 tolerance

// global variables
double kappa_cgs;                // kappa from Cook et al. 1994
double P_B_surf_cgs;             // baryonic pressure at surface




// declaring function to change precision of to_string
string to_string_with_precision(const double& value);

// declaring functions for sorting data
void triad_zip(const vector<double>& X, const vector<double>& Y, const vector<double>& Z, vector<tuple<double, double, double>>& zipped);
void triad_unzip(const vector<tuple<double, double, double>>& zipped, vector<double>& X, vector<double>& Y, vector<double>& Z);
void triad_sort(vector<double>& X, vector<double>& Y, vector<double>& Z);
void quartet_zip(const vector<double>& W, const vector<double>& X, const vector<double>& Y, const vector<double>& Z, vector<tuple<double, double, double, double>>& zipped);
void quartet_unzip(const vector<tuple<double, double, double, double>>& zipped, vector<double>& W, vector<double>& X, vector<double>& Y, vector<double>& Z);
void quartet_sort(vector<double>& W, vector<double>& X, vector<double>& Y, vector<double>& Z);

// declaring functions to calculate h_B
void h_B_integrate_trapezoid(vector<double>& epsilon_B_EOSfile_cgs, vector<double>& P_B_EOSfile_cgs, vector<double>& h_B_cgs);
double h_B_integrand(vector<double>& log10epsilon_B_EOSfile_cgs, vector<double>& log10P_B_EOSfile_cgs, const double& log10P_B_cgs, int *n_nearest_pt);
void h_B_integrate(vector<double>& log10epsilon_B_EOSfile_cgs, vector<double>& log10P_B_EOSfile_cgs, vector<double>& log10h_B_EOSfile_cgs, int *n_nearest_pt);

// declaring functions for root finding accuracy
bool P_B_c_rootfind_termination(const double& rootmin, const double& rootmax);
bool x_rootfind_termination(const double& rootmin, const double& rootmax);

// declaring functions for interpolating epsilon, P, h and n
//tk::spline log10h_B_cgs_of_log10epsilon_B_cgs;
//tk::spline log10P_B_cgs_of_log10epsilon_B_cgs;
//tk::spline log10epsilon_B_cgs_of_log10h_B_cgs;
//tk::spline log10epsilon_B_of_log10h_B;
double epsilon_B(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& h_B, int *n_nearest_pt);
//tk::spline log10P_B_cgs_of_log10h_B_cgs;
double P_B(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& h_B, int *n_nearest_pt);
//tk::spline log10n_B_cgs_of_log10h_B_cgs;
double n_B(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& h_B, int *n_nearest_pt);
//tk::spline log10P_B_cgs_of_log10epsilon_B_cgs;
double epsilon_D_EOS(const double& m_D, const double& x, const double& y, const string& units);
double P_D_EOS(const double& m_D, const double& x, const double& y, const string& units);
double h_D_EOS(const double& m_D, const double& x, const double& y, const string& units);
double epsilon_D(const double& h_D, const double& m_D, const double& y, const string& units);
double P_D(const double& h_D, const double& m_D, const double& y, const string& units);
double n_D(const double& h_D, const double& m_D, const double& y);

// declaring functions for the differential equations to be solved
double dr2dPhi(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y);
double dr2dh(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt);
double dM_Bdh(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt);
double dM_Ddh(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt);

// declaring function for finding dr_max
//double dr_max_region(const vector<double>& r);

// declaring function for RK4 method
void RK4_TOV(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, vector<double>& r2, vector<double>& h_B, vector<double>& M_B, vector<double>& h_D, vector<double>& M_D, const double& h_B_c, const double& h_D_c, const double& m_D, const double& y);
// declaring functions for RKF45 method
//void RKF45_TOV(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10P_B_EOSfile, vector<double>& r, vector<double>& P_B, vector<double>& M_B, vector<double>& P_D, vector<double>& M_D, const double& m_D, const double& y, const double& epsilon_B_c, const double& epsilon_D_c);

// declaring function to integrate Phi using trapezoid rule
double dPhidr(const double& r, const double& P_B, const double& M_B, const double& P_D, const double& M_D, const double& epsilon_B_c, const double& P_B_c, const double& epsilon_D_c, const double& P_D_c);
void metric_integrate_trapezoid(vector<double>& r, vector<double>& P_B, vector<double>& M_B, vector<double>& P_D, vector<double>& M_D, vector<double>& Phi, vector<double>& g, vector<double>& f, const double& epsilon_B_c, const double& epsilon_D_c, const int& R_B_index, const int& R_D_index);

//// declaring function to integrate N using trapezoid rule
//double dN_Bdr(const vector<double>& log10P_B_EOSfile, const vector<double>& log10n_B_EOSfile, const double& r, const double& P_B, const double& M_B, const double& M_D, const double& epsilon_B_c, const double& n_B_c, const double& epsilon_D_c);
//void N_B_integrate_trapezoid(const vector<double>& log10P_B_EOSfile, const vector<double>& log10n_B_EOSfile, vector<double>& r, vector<double>& P_B, vector<double>& M_B, vector<double>& N_B, vector<double>& M_D, const double& epsilon_B_c, const double& n_B_c, const double& epsilon_D_c, const int& R_B_index);
//double dN_Ddr(const double& r, const double& M_B, const double& P_D, const double& M_D, const double& m_D, const double& y, const double& epsilon_B_c, const double& epsilon_D_c, const double& n_D_c);
//void N_D_integrate_trapezoid(vector<double>& r, vector<double>& M_B, vector<double>& P_D, vector<double>& M_D, vector<double>& N_D, const double& m_D, const double& y, const double& epsilon_B_c, const double& epsilon_D_c, const double& n_D_c, const int& R_D_index);

// defining function to produce a single DANS solution
//vector<double> DANS_generator(const double& epsilon_B_c_cgs, const double& epsilon_D_c_cgs, const double& m_D_cgs, const double& y);
void DANS_generator(const double& epsilon_B_c_cgs, const double& epsilon_D_c_cgs, const double& m_D_cgs, const double& y);

//void hunt(const vector<double>&  xx, int n, double x, int *jlo);
void hunt(const vector<double>& xx, int n, double x, int *jlo);

//double interp(const vector<double>&  xp,
//              const vector<double>&  yp,
//              int    np,
//              double xb,
//              int n_nearest_pt);
double interp(const vector<double>& xp,
              const vector<double>& yp,
              int    np ,
              double xb,
              int    *n_nearest_pt);

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))




// file names
const string B_EOS_filename = "/Users/Tanjih/Documents/MSc-PhD_Thesis/Codes/eos-master/eos" + B_EOS;          // baryonic EOS
const string DANS_data_filename = "/Users/Tanjih/Documents/MSc-PhD_Thesis/Codes/non-rotating_DANS-selflensing_v3/DANS_enthalpy_RK4_interp_" + B_EOS + "_" + to_string_with_precision(m_D_MeVfm) + "MeV" + to_string_with_precision(y) + "_" + to_string_with_precision(epsilon_B_c_MeVfm) + "_" + to_string_with_precision(epsilon_D_c_MeVfm) + ".txt"; // DANS
const string enthalpy_data_filename = "/Users/Tanjih/Documents/MSc-PhD_Thesis/Codes/non-rotating_DANS_v6/enthalpy.txt"; // DANS
const string x_data_filename = "/Users/Tanjih/Documents/MSc-PhD_Thesis/Codes/non-rotating_DANS-selflensing_v3/x_enthalpy_RK4_" + B_EOS + "_" + to_string_with_precision(m_D_MeVfm) + "MeV" + to_string_with_precision(y) + ".txt"; // DANS




// cgs EOS file data
vector<double> log10epsilon_B_EOSfile_cgs;
vector<double> log10P_B_EOSfile_cgs;
vector<double> log10h_B_EOSfile_cgs;
vector<double> log10n_B_EOSfile_cgs;

int main()
{
    /* Time function returns the time since the
     Epoch(jan 1 1970). Returned time is in seconds. */
    time_t start, end;
     
    /* You can call it like this : start = time(NULL);
     in both the way start contain total time in seconds
     since the Epoch. */
    time(&start);
     
    // unsync the I/O of C and C++.
    ios_base::sync_with_stdio(false);
    
    cout << setprecision(12);

    
    // read data from EOS file and save log10 values
    ifstream B_EOS_file;
    B_EOS_file.open(B_EOS_filename);

    if(B_EOS_file.is_open())
    {
        if(B_EOS_CompOSE_flag == 0)
        {
            string nlines, epsilon_B_cgs, P_B_cgs, h_B_cgs, n_B_cgs; // rns B_EOS_file column headers
            B_EOS_file >> nlines;
            while(B_EOS_file >> epsilon_B_cgs >> P_B_cgs >> h_B_cgs >> n_B_cgs)
            {
                // take log10 and store
                log10epsilon_B_EOSfile_cgs.push_back(log10(stod(epsilon_B_cgs)));
                log10P_B_EOSfile_cgs.push_back(log10(stod(P_B_cgs)));
                log10h_B_EOSfile_cgs.push_back(log10(stod(h_B_cgs)));
                log10n_B_EOSfile_cgs.push_back(log10(stod(n_B_cgs)));
            }
        } else if(B_EOS_CompOSE_flag == 1)
        {
            string T_MeVfm, n_B_MeVfm, Y_q, epsilon_B_MeVfm, P_B_MeVfm, H_B_MeVfm; // CompOSE B_EOS_file column headers
            while(B_EOS_file >> T_MeVfm >> n_B_MeVfm >> Y_q >> epsilon_B_MeVfm >> P_B_MeVfm)
            {
                // convert from MeVfm to cgs, take log10 and store
                log10epsilon_B_EOSfile_cgs.push_back(log10(stod(epsilon_B_MeVfm)*MeV_cgs/(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0))));
                log10P_B_EOSfile_cgs.push_back(log10(stod(P_B_MeVfm)*MeV_cgs/pow(fm_cgs, 3.0)));
                log10n_B_EOSfile_cgs.push_back(log10(stod(n_B_MeVfm)/(pow(fm_cgs, 3.0))));
            }
            
            triad_sort(log10P_B_EOSfile_cgs, log10epsilon_B_EOSfile_cgs, log10n_B_EOSfile_cgs); // sort baryonic pressure in increasing order
            
            log10h_B_EOSfile_cgs.push_back(log10(h_B_surf_cgs));
            int n_nearest = log10P_B_EOSfile_cgs.size()/2;
            h_B_integrate(log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10h_B_EOSfile_cgs, &n_nearest);
        } else
        {
            cout << "Invalid B_EOS_CompOSE_flag" << endl;
        }
    }
    B_EOS_file.close();
    
    P_B_surf_cgs = pow(10.0, log10P_B_EOSfile_cgs[0]); // baryonic pressure at surface


    DANS_generator(epsilon_B_c_cgs, epsilon_D_c_cgs, m_D_cgs, y);
//    DANS_generator(1.8582e15, epsilon_D_c_cgs, m_D_cgs, y);



    // Recording end time.
    time(&end);
     
    // Calculating total time taken by the program.
    double time_taken = double(end - start);
    cout << "Time taken by program is : " << fixed
        << time_taken << setprecision(5);
    cout << " sec " << endl;



    return 0;
}




// defining function to change precision of to_string
string to_string_with_precision(const double& value)
{
    ostringstream out;
    int decimal_pts = 0;
    int count = 0;
    
    while(count >= 0)
    {
        ostringstream out;
        out.precision(decimal_pts);
        out << scientific << value;
        if(stod(out.str()) == value)
        {
            break;
        } else
        {
            decimal_pts += 1;
            count += 1;
        }
    }
    
    out.precision(decimal_pts);
    out << scientific << value;
    
    return out.str();
}




// declaring functions for sorting data
void triad_zip(const vector<double>& X, const vector<double>& Y, const vector<double>& Z, vector<tuple<double, double, double>>& zipped)
{
    for(int i = 0; i < X.size(); i++)
    {
        zipped.push_back(make_tuple(X[i], Y[i], Z[i]));
    }
}

void triad_unzip(const vector<tuple<double, double, double>>& zipped, vector<double>& X, vector<double>& Y, vector<double>& Z)
{
    for(int i = 0; i < X.size(); i++)
    {
        X[i] = get<0>(zipped[i]);
        Y[i] = get<1>(zipped[i]);
        Z[i] = get<2>(zipped[i]);
    }
}

void triad_sort(vector<double>& X, vector<double>& Y, vector<double>& Z)
{
    // Zip the vectors together
    vector<tuple<double, double, double>> zipped;
    triad_zip(X, Y, Z, zipped);

    // Sort the vector of triads
    sort(begin(zipped), end(zipped),
         [&](const tuple<double, double, double>& a, const tuple<double, double, double>& b)
         {
             return get<0>(a) < get<0>(b);
         });

    // Write the sorted quartets back to the original vectors
    triad_unzip(zipped, X, Y, Z);
}


void quartet_zip(const vector<double>& W, const vector<double>& X, const vector<double>& Y, const vector<double>& Z, vector<tuple<double, double, double, double>>& zipped)
{
    for(int i = 0; i < W.size(); i++)
    {
        zipped.push_back(make_tuple(W[i], X[i], Y[i], Z[i]));
    }
}

void quartet_unzip(const vector<tuple<double, double, double, double>>& zipped, vector<double>& W, vector<double>& X, vector<double>& Y, vector<double>& Z)
{
    for(int i = 0; i < X.size(); i++)
    {
        W[i] = get<0>(zipped[i]);
        X[i] = get<1>(zipped[i]);
        Y[i] = get<2>(zipped[i]);
        Z[i] = get<3>(zipped[i]);
    }
}

void quartet_sort(vector<double>& W, vector<double>& X, vector<double>& Y, vector<double>& Z)
{
    // Zip the vectors together
    vector<tuple<double, double, double, double>> zipped;
    quartet_zip(W, X, Y, Z, zipped);

    // Sort the vector of quartets
    sort(begin(zipped), end(zipped),
         [&](const tuple<double, double, double, double>& a, const tuple<double, double, double, double>& b)
         {
             return get<0>(a) < get<0>(b);
         });

    // Write the sorted quartets back to the original vectors
    quartet_unzip(zipped, W, X, Y, Z);
}





// defining function to calculate h_B
void h_B_integrate_trapezoid(vector<double>& epsilon_B_EOSfile_cgs, vector<double>& P_B_EOSfile_cgs, vector<double>& h_B_cgs)
{
    double h_B_integral = 0.0;
    for(int i = 1; i < epsilon_B_EOSfile_cgs.size(); i++)
    {
        h_B_integral += (pow(c_cgs, 2.0)/(epsilon_B_EOSfile_cgs[i - 1]*pow(c_cgs, 2.0) + P_B_EOSfile_cgs[i - 1]) + pow(c_cgs, 2.0)/(epsilon_B_EOSfile_cgs[i]*pow(c_cgs, 2.0) + P_B_EOSfile_cgs[i]))*(P_B_EOSfile_cgs[i] - P_B_EOSfile_cgs[i - 1])/2.0;
        h_B_cgs.push_back(h_B_integral);
    }
}

double h_B_integrand(vector<double>& log10epsilon_B_EOSfile_cgs, vector<double>& log10P_B_EOSfile_cgs, const double& log10P_B_cgs, int *n_nearest_pt)
{
    return (1.0/(1.0 + pow(10.0, interp(log10P_B_EOSfile_cgs, log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs.size(), log10P_B_cgs, n_nearest_pt) + log10(pow(c_cgs, 2.0)) - log10P_B_cgs)));
}

void h_B_integrate(vector<double>& log10epsilon_B_EOSfile_cgs, vector<double>& log10P_B_EOSfile_cgs, vector<double>& log10h_B_EOSfile_cgs, int *n_nearest_pt)
{
    int n = 16003;
    double log10P_B_cgs[n];
    double log10P_B_cgs_val;
    double dlog10P_B_cgs;
    double h_B_integral;
    for(int i = 1; i < log10epsilon_B_EOSfile_cgs.size(); i++)
    {
        for(int j = 0; j < n; j++)
        {
            log10P_B_cgs_val = log10P_B_EOSfile_cgs[0] + (log10P_B_EOSfile_cgs[i] - log10P_B_EOSfile_cgs[0])*j/(n - 1);
            if(log10P_B_cgs_val < log10P_B_EOSfile_cgs[0])
            {
                log10P_B_cgs[j] = log10P_B_EOSfile_cgs[0];
            } else if(log10P_B_cgs_val > log10P_B_EOSfile_cgs.back())
            {
                log10P_B_cgs[j] = log10P_B_EOSfile_cgs.back();
            } else
            {
                log10P_B_cgs[j] = log10P_B_cgs_val;
            }
        }
        
        dlog10P_B_cgs = (log10P_B_EOSfile_cgs[i] - log10P_B_EOSfile_cgs[0])/(n - 1);
        
        h_B_integral = 0.0;
        for(int j = 0; j < n - 2; j += 2)
        {
            h_B_integral += (dlog10P_B_cgs/3.0)*(h_B_integrand(log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10P_B_cgs[j], n_nearest_pt) + 4.0*h_B_integrand(log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10P_B_cgs[j + 1], n_nearest_pt) + h_B_integrand(log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10P_B_cgs[j + 2], n_nearest_pt));
        }
        
        log10h_B_EOSfile_cgs.push_back(log10(pow(c_cgs, 2.0)*log(10.0)*h_B_integral));
    }
}




// defining functions for root finding
bool P_B_c_rootfind_termination(const double& rootmin, const double& rootmax)
{
  return abs((rootmax - rootmin)/rootmin) <= P_B_c_rootfind_accuracy;
}

bool x_rootfind_termination(const double& rootmin, const double& rootmax)
{
  return abs((rootmax - rootmin)/rootmin) <= x_rootfind_accuracy;
}




// defining functions for epsilon, P and n
double epsilon_B(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& h_B, int *n_nearest_pt)
{
    if(h_B >= pow(10.0, log10h_B_EOSfile[0]))
    {
        return pow(10.0, interp(log10h_B_EOSfile, log10epsilon_B_EOSfile, log10h_B_EOSfile.size(), log10(h_B), n_nearest_pt));
    } else
    {
        return epsilon_B_surf_zero;
    }
}

double P_B(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& h_B, int *n_nearest_pt)
{
    if(h_B >= pow(10.0, log10h_B_EOSfile[0]))
    {
        return pow(10.0, interp(log10h_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile.size(), log10(h_B), n_nearest_pt));
    } else
    {
        return P_B_surf_zero;
    }
}

double n_B(const vector<double>& log10h_B_EOSfile, const vector<double>& log10n_B_EOSfile, const double& h_B, int *n_nearest_pt)
{
    if(h_B >= pow(10.0, log10h_B_EOSfile[0]))
    {
        return pow(10.0, interp(log10h_B_EOSfile, log10n_B_EOSfile, log10h_B_EOSfile.size(), log10(h_B), n_nearest_pt));
    } else
    {
        return n_B_surf_zero;
    }
}

double epsilon_D_EOS(const double& m_D, const double& x, const double& y, const string& units)
{
    if(units == "dimensionless")
    {
        return (pow(m_D, 4.0)/(pow(hbar, 3.0)*pow(M_PI, 2.0)))*((1.0/8.0)*((2.0*pow(x, 3.0) + x)*sqrt(1.0 + pow(x, 2.0)) - asinh(x)) + pow(y, 2.0)*pow(x, 6.0)/(9.0*pow(M_PI, 2.0)));
    } else if(units == "cgs")
    {
        return (pow(m_D, 4.0)*pow(c_cgs/hbar_cgs, 3.0)/pow(M_PI, 2.0))*((1.0/8.0)*((2.0*pow(x, 3.0) + x)*sqrt(1.0 + pow(x, 2.0)) - asinh(x)) + pow(y, 2.0)*pow(x, 6.0)/(9.0*pow(M_PI, 2.0)));
    }
    else
    {
        cout << "Error" << endl;
        return -9999;
    }
}

double P_D_EOS(const double& m_D, const double& x, const double& y, const string& units)
{
    if(units == "dimensionless")
    {
        return (pow(m_D, 4.0)/(3.0*pow(hbar, 3.0)*pow(M_PI, 2.0)))*((1.0/8.0)*((2.0*pow(x, 3.0) - 3.0*x)*sqrt(1.0 + pow(x, 2.0)) + 3.0*asinh(x)) + pow(y, 2.0)*pow(x, 6.0)/(3.0*pow(M_PI, 2.0)));
    } else if(units == "cgs")
    {
        return (pow(m_D, 4.0)*pow(c_cgs, 5.0)/(3.0*pow(hbar_cgs, 3.0)*pow(M_PI, 2.0)))*((1.0/8.0)*((2.0*pow(x, 3.0) - 3.0*x)*sqrt(1.0 + pow(x, 2.0)) + 3.0*asinh(x)) + pow(y, 2.0)*pow(x, 6.0)/(3.0*pow(M_PI, 2.0)));
    }
    else
    {
        cout << "Error" << endl;
        return -9999;
    }
}

double h_D_EOS(const double& m_D, const double& x, const double& y, const string& units)
{
    if(units == "dimensionless")
    {
        return log(sqrt(1.0 + pow(x, 2.0)) + 2.0*pow(x, 3.0)*pow(y, 2.0)/(3.0*pow(M_PI, 2.0)));
    } else if(units == "cgs")
    {
        return pow(c_cgs, 2.0)*log(sqrt(1.0 + pow(x, 2.0)) + 2.0*pow(x, 3.0)*pow(y, 2.0)/(3.0*pow(M_PI, 2.0)));
    }
    else
    {
        cout << "Error" << endl;
        return -9999;
    }
}

double epsilon_D(const double& h_D, const double& m_D, const double& y, const string& units)
{
    if(h_D > h_D_surf_zero)
    {
        pair<double, double> x_interval = boost::math::tools::bisect([m_D, y, h_D, units](double x_rootfind){return h_D_EOS(m_D, x_rootfind, y, units) - h_D;},
                                                                     x_min, x_max, x_rootfind_termination);
        double x = (x_interval.first + x_interval.second)/2.0;
        //const double k_F = pow(3.0*pow(M_PI, 2.0)*n_D, 1.0/3.0) // Fermi momentum
        //const double x = k_F/m_D                                // x parameter from Miao et al. 2022
        //const double y = m_D/m_I                                // y parameter from Miao et al. 2022
        return epsilon_D_EOS(m_D, x, y, units);
    } else
    {
        return epsilon_D_surf_zero;
    }
}

double P_D(const double& h_D, const double& m_D, const double& y, const string& units)
{
    if(h_D > h_D_surf_zero)
    {
        pair<double, double> x_interval = boost::math::tools::bisect([m_D, y, h_D, units](double x_rootfind){return h_D_EOS(m_D, x_rootfind, y, units) - h_D;},
                                                                     x_min, x_max, x_rootfind_termination);
        double x = (x_interval.first + x_interval.second)/2.0;
        //const double k_F = pow(3.0*pow(M_PI, 2.0)*n_D, 1.0/3.0) // Fermi momentum
        //const double x = k_F/m_D                                // x parameter from Miao et al. 2022
        //const double y = m_D/m_I                                // y parameter from Miao et al. 2022
        return P_D_EOS(m_D, x, y, units);
    } else
    {
        return P_D_surf_zero;
    }
}

//double n_D(const double& h_D, const double& m_D, const double& y, const string& units)
//{
//    if(h_D > h_D_surf_zero)
//    {
//        pair<double, double> x_interval = boost::math::tools::bisect([m_D, y, h_D, units](double x_rootfind){return h_D_EOS(m_D, x_rootfind, y, units) - h_D;}, x_min, x_max, x_rootfind_termination);
//        double x = (x_interval.first + x_interval.second)/2.0;
//        //const double k_F = pow(3.0*pow(M_PI, 2.0)*n_D, 1.0/3.0) // Fermi momentum
//        //const double x = k_F/m_D                                // x parameter from Miao et al. 2022
//        //const double y = m_D/m_I                                // y parameter from Miao et al. 2022
//        return pow(m_D*x/hbar, 3.0)/(3.0*pow(M_PI, 2.0));
//    } else
//    {
//        return n_D_surf;
//    }
//}




// defining functions of the differential equations to be solved
double dr2dPhi(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt)
{
    return 2.0*r2*(sqrt(r2) - 2.0*(M_B + M_D))/(M_B + M_D + 4.0*M_PI*pow(r2, 1.5)*(P_B(log10P_B_EOSfile, log10h_B_EOSfile, h_B, n_nearest_pt) + P_D(h_D, m_D, y, "dimensionless")));
}

double dr2dh(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt)
{
    return -dr2dPhi(log10P_B_EOSfile, log10h_B_EOSfile, r2, h_B, M_B, h_D, M_D, m_D, y, n_nearest_pt);
}

double dM_Bdh(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt)
{
    return 2.0*M_PI*sqrt(r2)*epsilon_B(log10epsilon_B_EOSfile, log10h_B_EOSfile, h_B, n_nearest_pt)*dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, r2, h_B, M_B, h_D, M_D, m_D, y, n_nearest_pt);
}

//double dN_Bdr(const double& r, const double& h_B, const double& M_B, const double& M_D, const double& epsilon_B_c, const double& n_B_c, const double& epsilon_D_c)
//{
//    if(r > r_powerseries_max)
//    {
//        return 4.0*M_PI*pow(r, 2.0)*n_B(h_B)/sqrt(1.0 - 2.0*(M_B + M_D)/r);
//    } else
//    {
//        return 4.0*M_PI*n_B_c*pow(r, 2.0)/sqrt(1.0 - 8.0*M_PI*(epsilon_B_c + epsilon_D_c)*pow(r, 2.0)/3.0);
//    }
//}

double dM_Ddh(const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, const double& r2, const double& h_B, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, int *n_nearest_pt)
{
    return 2.0*M_PI*sqrt(r2)*epsilon_D(h_D, m_D, y, "dimensionless")*dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, r2, h_B, M_B, h_D, M_D, m_D, y, n_nearest_pt);
}

//double dN_Ddr(const double& r, const double& M_B, const double& h_D, const double& M_D, const double& m_D, const double& y, const double& epsilon_B_c, const double& epsilon_D_c, const double& n_D_c)
//{
//    if(r > r_powerseries_max)
//    {
//        return 4.0*M_PI*pow(r, 2.0)*n_D(h_D, m_D, y)/sqrt(1.0 - 2.0*(M_B + M_D)/r);
//    } else
//    {
//        return 4.0*M_PI*n_D_c*pow(r, 2.0)/sqrt(1.0 - 8.0*M_PI*(epsilon_B_c + epsilon_D_c)*pow(r, 2.0)/3.0);
//    }
//}




// defining function for RK4 method
void RK4_TOV(const vector<double>& log10epsilon_B_EOSfile, const vector<double>& log10P_B_EOSfile, const vector<double>& log10h_B_EOSfile, vector<double>& r2, vector<double>& h_B, vector<double>& M_B, vector<double>& h_D, vector<double>& M_D, const double& h_B_c, const double& h_D_c, const double& m_D, const double& y)
{
    double k1_r2input, k1_r2, k2_r2input, k2_r2, k3_r2input, k3_r2, k4_r2input, k4_r2, r2input;
    double k1_h_Binput, k2_h_Binput, k3_h_Binput, k4_h_Binput, h_Binput;
    double k1_M_Binput, k1_M_B, k2_M_Binput, k2_M_B, k3_M_Binput, k3_M_B, k4_M_Binput, k4_M_B, M_Binput;
    double k1_h_Dinput, k2_h_Dinput, k3_h_Dinput, k4_h_Dinput, h_Dinput;
    double k1_M_Dinput, k1_M_D, k2_M_Dinput, k2_M_D, k3_M_Dinput, k3_M_D, k4_M_Dinput, k4_M_D, M_Dinput;
    
    int n_nearest = log10h_B_EOSfile.size()/2;

    while(h_B.back() >= h_B_surf || h_D.back() > h_D_surf_zero)
    {
//        cout << "r = " << sqrt(r2.back())*sqrt(kappa_cgs)/km_cgs << " km" << endl;
//        cout << h_B_vec.back()*pow(c_cgs, 2.0) << endl;

        if(r2.back() > r_powerseries_max*r_powerseries_max)
        {
            k1_r2input = r2.back();
            k1_h_Binput = h_B.back();
            k1_M_Binput = M_B.back();
            k1_h_Dinput = h_D.back();
            k1_M_Dinput = M_D.back();
            if(k1_h_Binput < h_B_surf)
            {
                k1_h_Binput = h_B_surf_zero;
                k1_M_Binput = M_B.back();
            }
            if(k1_h_Dinput <= h_D_surf_zero)
            {
                k1_h_Dinput = h_D_surf_zero;
                k1_M_Dinput = M_D.back();
            }
            k1_r2 = dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, k1_r2input, k1_h_Binput, k1_M_Binput, k1_h_Dinput, k1_M_Dinput, m_D, y, &n_nearest);
            k1_M_B = dM_Bdh(log10epsilon_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile, k1_r2input, k1_h_Binput, k1_M_Binput, k1_h_Dinput, k1_M_Dinput, m_D, y, &n_nearest);
            k1_M_D = dM_Ddh(log10P_B_EOSfile, log10h_B_EOSfile, k1_r2input, k1_h_Binput, k1_M_Binput, k1_h_Dinput, k1_M_Dinput, m_D, y, &n_nearest);
            //cout << k1_r << endl;

            k2_r2input = r2.back() + k1_r2*dh/2.0;
            k2_h_Binput = h_B.back() + dh/2.0;
            k2_M_Binput = M_B.back() + k1_M_B*dh/2.0;
            k2_h_Dinput = h_D.back() + dh/2.0;
            k2_M_Dinput = M_D.back() + k1_M_D*dh/2.0;
            if(k2_h_Binput < h_B_surf)
            {
                k2_h_Binput = h_B_surf_zero;
                k2_M_Binput = M_B.back();
            }
            if(k2_h_Dinput <= h_D_surf_zero)
            {
                k2_h_Dinput = h_D_surf_zero;
                k2_M_Dinput = M_D.back();
            }
            k2_r2 = dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, k2_r2input, k2_h_Binput, k2_M_Binput, k2_h_Dinput, k2_M_Dinput, m_D, y, &n_nearest);
            k2_M_B = dM_Bdh(log10epsilon_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile, k2_r2input, k2_h_Binput, k2_M_Binput, k2_h_Dinput, k2_M_Dinput, m_D, y, &n_nearest);
            k2_M_D = dM_Ddh(log10P_B_EOSfile, log10h_B_EOSfile, k2_r2input, k2_h_Binput, k2_M_Binput, k2_h_Dinput, k2_M_Dinput, m_D, y, &n_nearest);


            k3_r2input = r2.back() + k2_r2*dh/2.0;
            k3_h_Binput = h_B.back() + dh/2.0;
            k3_M_Binput = M_B.back() + k2_M_B*dh/2.0;
            k3_h_Dinput = h_D.back() + dh/2.0;
            k3_M_Dinput = M_D.back() + k2_M_D*dh/2.0;
            if(k3_h_Binput < h_B_surf)
            {
                k3_h_Binput = h_B_surf_zero;
                k3_M_Binput = M_B.back();
            }
            if(k3_h_Dinput <= h_D_surf_zero)
            {
                k3_h_Dinput = h_D_surf_zero;
                k3_M_Dinput = M_D.back();
            }
            k3_r2 = dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, k3_r2input, k3_h_Binput, k3_M_Binput, k3_h_Dinput, k3_M_Dinput, m_D, y, &n_nearest);
            k3_M_B = dM_Bdh(log10epsilon_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile, k3_r2input, k3_h_Binput, k3_M_Binput, k3_h_Dinput, k3_M_Dinput, m_D, y, &n_nearest);
            k3_M_D = dM_Ddh(log10P_B_EOSfile, log10h_B_EOSfile, k3_r2input, k3_h_Binput, k3_M_Binput, k3_h_Dinput, k3_M_Dinput, m_D, y, &n_nearest);


            k4_r2input = r2.back() + k3_r2*dh;
            k4_h_Binput = h_B.back() + dh;
            k4_M_Binput = M_B.back() + k3_M_B*dh;
            k4_h_Dinput = h_D.back() + dh;
            k4_M_Dinput = M_D.back() + k2_M_D*dh;
            if(k4_h_Binput < h_B_surf)
            {
                k4_h_Binput = h_B_surf_zero;
                k4_M_Binput = M_B.back();
            }
            if(k4_h_Dinput <= h_D_surf_zero)
            {
                k4_h_Dinput = h_D_surf_zero;
                k4_M_Dinput = M_D.back();
            }
            k4_r2 = dr2dh(log10P_B_EOSfile, log10h_B_EOSfile, k4_r2input, k4_h_Binput, k4_M_Binput, k4_h_Dinput, k4_M_Dinput, m_D, y, &n_nearest);
            k4_M_B = dM_Bdh(log10epsilon_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile, k4_r2input, k4_h_Binput, k4_M_Binput, k4_h_Dinput, k4_M_Dinput, m_D, y, &n_nearest);
            k4_M_D = dM_Ddh(log10P_B_EOSfile, log10h_B_EOSfile, k4_r2input, k4_h_Binput, k4_M_Binput, k4_h_Dinput, k4_M_Dinput, m_D, y, &n_nearest);


            r2input = r2.back() + (k1_r2 + 2.0*k2_r2 + 2.0*k3_r2 + k4_r2)*dh/6.0;
            h_Binput = h_B.back() + dh;
            M_Binput = M_B.back() + (k1_M_B + 2.0*k2_M_B + 2.0*k3_M_B + k4_M_B)*dh/6.0;
            h_Dinput = h_D.back() + dh;
            M_Dinput = M_D.back() + (k1_M_D + 2.0*k2_M_D + 2.0*k3_M_D + k4_M_D)*dh/6.0;
            if(k1_h_Binput < h_B_surf || k2_h_Binput < h_B_surf || k3_h_Binput < h_B_surf || k4_h_Binput < h_B_surf || h_Binput < h_B_surf)
            {
                h_Binput = h_B_surf_zero;
                M_Binput = M_B.back();
            }
            if(k1_h_Dinput <= h_D_surf_zero || k2_h_Dinput <= h_D_surf_zero || k3_h_Dinput <= h_D_surf_zero || k4_h_Dinput <= h_D_surf_zero || h_Dinput <= h_D_surf_zero)
            {
                h_Dinput = h_D_surf_zero;
                M_Dinput = M_D.back();
            }
        } else
        {
            h_Binput = h_B.back() + dh;
            h_Dinput = h_D.back() + dh;
            r2input = 3.0*(h_B_c + h_D_c - h_Binput - h_Dinput)/(2.0*M_PI*(3.0*(P_B(log10P_B_EOSfile, log10h_B_EOSfile, h_B_c, &n_nearest) + P_D(h_D_c, m_D, y, "dimensionless")) + epsilon_B(log10epsilon_B_EOSfile, log10h_B_EOSfile, h_B_c, &n_nearest) + epsilon_D(h_D_c, m_D, y, "dimensionless")));
            M_Binput = 4.0*M_PI*epsilon_B(log10epsilon_B_EOSfile, log10h_B_EOSfile, h_B_c, &n_nearest)*pow(r2input, 1.5)/3.0;
            M_Dinput = 4.0*M_PI*epsilon_D(h_D_c, m_D, y, "dimensionless")*pow(r2input, 1.5)/3.0;

            if(h_Binput < h_B_surf)
            {
                h_Binput = h_B_surf_zero;
                M_Binput = M_B.back();
            }
            if(h_Dinput <= h_D_surf_zero)
            {
                h_Dinput = h_D_surf_zero;
                M_Dinput = M_D.back();
            }
        }


        r2.push_back(r2input);
        h_B.push_back(h_Binput);
        M_B.push_back(M_Binput);
        h_D.push_back(h_Dinput);
        M_D.push_back(M_Dinput);
    }
}




// defining function to integrate Phi using trapezoid rule
double dPhidr(const double& r, const double& P_B, const double& M_B, const double& P_D, const double& M_D, const double& epsilon_B_c, const double& P_B_c, const double& epsilon_D_c, const double& P_D_c)
{
    if(r < r_powerseries_max)
    {
        return 4.0*M_PI*((epsilon_B_c + epsilon_D_c)/3.0 + P_B_c + P_D_c)*r/(1.0 - 8.0*M_PI*(epsilon_B_c + epsilon_D_c)*pow(r, 2.0)/3.0);
    } else
    {
        return (M_B + M_D + 4.0*M_PI*pow(r, 3.0)*(P_B + P_D))/(r*(r - 2.0*(M_B + M_D)));
    }
}

void metric_integrate_trapezoid(vector<double>& r, vector<double>& P_B, vector<double>& M_B, vector<double>& P_D, vector<double>& M_D, vector<double>& Phi, vector<double>& g, vector<double>& f, const double& epsilon_B_c, const double& epsilon_D_c, const int& R_B_index, const int& R_D_index)
{
    double Phi_integral = 0.0;
    for(int i = max(R_B_index, R_D_index) - 1; i >= 0; i--)
    {
        Phi_integral += (dPhidr(r[i], P_B[i], M_B[i], P_D[i], M_D[i], epsilon_B_c, P_B[0], epsilon_D_c, P_D[0]) + dPhidr(r[i + 1], P_B[i + 1], M_B[i + 1], P_D[i + 1], M_D[i + 1], epsilon_B_c, P_B[0], epsilon_D_c, P_D[0]))*(r[i + 1] - r[i])/2.0;
        Phi.insert(Phi.begin(), Phi.back() - Phi_integral);
        g.insert(g.begin(), exp(2.0*Phi[0]));
        if(i != 0)
        {
            f.insert(f.begin(), pow(1.0 - 2.0*(M_B[i] + M_D[i])/r[i], -1.0));
        } else
        {
            f.insert(f.begin(), f_c);
        }
    }
}




//// defining function to integrate N using trapezoid rule
//void N_B_integrate_trapezoid(const vector<double>& log10P_B_EOSfile, const vector<double>& log10n_B_EOSfile, vector<double>& r, vector<double>& P_B, vector<double>& M_B, vector<double>& N_B, vector<double>& M_D, const double& epsilon_B_c, const double& n_B_c, const double& epsilon_D_c, const int& R_B_index)
//{
//    double N_B_integral = 0.0;
//    for(int i = 1; i <= R_B_index; i++)
//    {
//        N_B_integral += (dN_Bdr(log10P_B_EOSfile, log10n_B_EOSfile, r[i - 1], P_B[i - 1], M_B[i - 1], M_D[i - 1], epsilon_B_c, n_B_c, epsilon_D_c) + dN_Bdr(log10P_B_EOSfile, log10n_B_EOSfile, r[i], P_B[i], M_B[i], M_D[i], epsilon_B_c, n_B_c, epsilon_D_c))*(r[i] - r[i - 1])/2.0;
//        N_B.push_back(N_B_integral);
//    }
//}
//
//void N_D_integrate_trapezoid(vector<double>& r, vector<double>& M_B, vector<double>& P_D, vector<double>& M_D, vector<double>& N_D, const double& m_D, const double& y, const double& epsilon_B_c, const double& epsilon_D_c, const double& n_D_c, const int& R_D_index)
//{
//    double N_D_integral = 0.0;
//    for(int i = 1; i <= R_D_index; i++)
//    {
//        N_D_integral += (dN_Ddr(r[i - 1], M_B[i - 1], P_D[i - 1], M_D[i - 1], m_D, y, epsilon_B_c, epsilon_D_c, n_D_c) + dN_Ddr(r[i], M_B[i], P_D[i], M_D[i], m_D, y, epsilon_B_c, epsilon_D_c, n_D_c))*(r[i] - r[i - 1])/2.0;
//        N_D.push_back(N_D_integral);
//    }
//}




// declaring function to produce a single DANS solution
//vector<double> DANS_generator(const double& epsilon_B_c_cgs, const double& epsilon_D_c_cgs, const double& m_D_cgs, const double& y)
void DANS_generator(const double& epsilon_B_c_cgs, const double& epsilon_D_c_cgs, const double& m_D_cgs, const double& y)
{
    // find the correct central baryonic specific enthalpy
    cout << "Finding correct h_B_c_cgs..." << endl;
    int n_nearest = log10h_B_EOSfile_cgs.size()/2;
    double h_B_c_cgs = pow(10.0, interp(log10epsilon_B_EOSfile_cgs, log10h_B_EOSfile_cgs, log10h_B_EOSfile_cgs.size(), log10(epsilon_B_c_cgs), &n_nearest));
    cout << "h_B_c_cgs found = " << h_B_c_cgs << endl;
    
    // find the correct central baryonic pressure
//    cout << "Finding correct P_B_c_cgs..." << endl;
    double P_B_c_cgs = pow(10.0, interp(log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10P_B_EOSfile_cgs.size(), log10(epsilon_B_c_cgs), &n_nearest));
    cout << "P_B_c_cgs found = " << P_B_c_cgs << endl;
    cout << "Correct epsilon_B_c_cgs = " << epsilon_B_c_cgs << endl;
    cout << "epsilon_B_c_cgs found = " << epsilon_B(log10epsilon_B_EOSfile_cgs, log10h_B_EOSfile_cgs, h_B_c_cgs, &n_nearest) << endl;
    
    quartet_sort(log10h_B_EOSfile_cgs, log10epsilon_B_EOSfile_cgs, log10P_B_EOSfile_cgs, log10n_B_EOSfile_cgs); // sort baryonic specific enthalpy in increasing order

    // find the correct central DM specific enthalpy
//    cout << "Finding correct P_D_c_cgs..." << endl;
    pair<double, double> x_c_interval = boost::math::tools::bisect([epsilon_D_c_cgs, m_D_cgs, y](double x_c_rootfind){return epsilon_D_EOS(m_D_cgs, x_c_rootfind, y, "cgs") - epsilon_D_c_cgs;}, x_min, x_max, x_rootfind_termination);
    double x_c = (x_c_interval.first + x_c_interval.second)/2.0;
    double h_D_c_cgs = h_D_EOS(m_D_cgs, x_c, y, "cgs");
    cout << "h_D_c_cgs found = " << h_D_c_cgs << endl;

    // find the correct central DM pressure
//    cout << "Finding correct P_D_c_cgs..." << endl;
    double P_D_c_cgs = P_D_EOS(m_D_cgs, x_c, y, "cgs");
    cout << "x_c found = " << x_c << endl;
    cout << "P_D_c_cgs found = " << P_D_c_cgs << endl;
    cout << "Correct epsilon_D_c_cgs = " << epsilon_D_c_cgs << endl;
    cout << "epsilon_D_c_cgs found = " << epsilon_D_EOS(m_D_cgs, x_c, y, "cgs") << endl;
    
    
    
    // find the correct central DM specific enthalpy
//    cout << "Finding correct P_D_c_cgs..." << endl;
    pair<double, double> x_surf_interval = boost::math::tools::bisect([m_D_cgs, y](double x_c_rootfind){return h_D_EOS(m_D_cgs, x_c_rootfind, y, "cgs") - pow(10.0, log10h_B_EOSfile_cgs[0]);}, x_min, x_max, x_rootfind_termination);
    double x_surf = (x_surf_interval.first + x_surf_interval.second)/2.0;
    double h_D_surf_cgs = h_D_EOS(m_D_cgs, x_surf, y, "cgs");
    cout << "Correct h_D_surf_cgs = " << pow(10.0, log10h_B_EOSfile_cgs[0]) << endl;
    cout << "h_D_surf_cgs found = " << h_D_surf_cgs << endl;

    // find the correct central DM pressure
//    cout << "Finding correct P_D_c_cgs..." << endl;
    cout << "x_surf found = " << x_surf << endl;
//    cout << "P_D_c_cgs found = " << P_D_c_cgs << endl;
//    cout << "Corerect epsilon_D_c_cgs = " << epsilon_D_c_cgs << endl;
//    cout << "epsilon_D_c_cgs found = " << epsilon_D_EOS(m_D_cgs, x_c, y, "cgs") << endl;


    // find appropriate value of epsilon_0_cgs from Cook et al. 1994
    double epsilon_0_cgs;
    if(epsilon_D_c_cgs != epsilon_D_c_zero && epsilon_B_c_cgs != epsilon_B_c_zero)
    {
        epsilon_0_cgs = min({P_B_c_cgs, P_D_c_cgs})/pow(c_cgs, 2.0);
    } else if(epsilon_D_c_cgs == epsilon_D_c_zero)
    {
        epsilon_0_cgs = P_B_c_cgs/pow(c_cgs, 2.0);
    } else
    {
        epsilon_0_cgs = P_D_c_cgs/pow(c_cgs, 2.0);
    }
//    cout << "epsilon_0_cgs = " << epsilon_0_cgs << endl;
    kappa_cgs = pow(c_cgs, 2.0)/(G_cgs*epsilon_0_cgs); // kappa from Cook et al. 1994
//    cout << "kappa_cgs = " << kappa_cgs << endl;


    // unit conversion to dimensionless
    hbar = hbar_cgs*G_cgs/(kappa_cgs*pow(c_cgs, 3.0));                       // reduced Planck's constant
    //M_T = M_T_cgs*G_cgs/(sqrt(kappa_cgs)*pow(c_cgs, 2.0));             // total mass of DANS
    //M_B = M_T_input*(1.0 - f_D_input);                                 // baryonic mass of DANS
    //M_D = M_T_input*f_D_input;                                         // DM mass of DANS
    //m_phi = G_cgs*m_phi_cgs/(pow(kappa_cgs, 1.0/2.0)*pow(c_cgs, 2.0)); // DM self-interaction mediator mass
    double m_D = m_D_cgs*G_cgs/(sqrt(kappa_cgs)*pow(c_cgs, 2.0));            // DM particle mass
    r_TOV_accurateregion_max = r_TOV_accurateregion_max_cgs/sqrt(kappa_cgs); // upper limit for high resolution region for solving TOV equations
    dr_TOV_min = dr_TOV_min_cgs/sqrt(kappa_cgs);                             // minimum grid spacing for solving TOV equations
//    dr_TOV_max = dr_TOV_max_cgs/sqrt(kappa_cgs);                             // maximum grid spacing for solving TOV equations
    r_powerseries_max = dr_TOV_min*r_powerseries_max_factor;                 // end point of power series
    dr_TOV_accurateregion_max = dr_TOV_min*dr_TOV_accurateregion_max_factor; // maximum grid spacing in high resolution region for solving TOV equations
    double epsilon_B_c = epsilon_B_c_cgs*kappa_cgs*G_cgs/pow(c_cgs, 2.0);    // baryonic central energy density
    double P_B_c = P_B_c_cgs*kappa_cgs*G_cgs/pow(c_cgs, 4.0);                // baryonic central pressure
    double h_B_c = h_B_c_cgs/pow(c_cgs, 2.0);                                // baryonic central specific enthalpy
    P_B_surf = P_B_surf_cgs*kappa_cgs*G_cgs/pow(c_cgs, 4.0);                 // baryonic pressure at surface
    double epsilon_D_c = epsilon_D_c_cgs*kappa_cgs*G_cgs/pow(c_cgs, 2.0);    // DM central energy densi
    double P_D_c = P_D_c_cgs*kappa_cgs*G_cgs/pow(c_cgs, 4.0);                // DM central pressure
    double h_D_c = h_D_c_cgs/pow(c_cgs, 2.0);                                // DM central specific enthalpy
    h_B_surf = h_B_surf_cgs/pow(c_cgs, 2.0);                          // baryonic specific enthalpy at surface

    if(epsilon_D_c_cgs != epsilon_D_c_zero && epsilon_B_c_cgs != epsilon_B_c_zero)
    {
        dh = -min(h_B_c, h_D_c)/1.0e4;
    } else if(epsilon_D_c_cgs == epsilon_D_c_zero)
    {
        dh = -h_B_c/1.0e4;
    } else
    {
        dh = -h_D_c/1.0e4;
    }

    cout << "epsilon_B_c_MeVfm = " << epsilon_B_c_cgs*(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0))/MeV_cgs << endl;
    cout << "epsilon_D_c_MeVfm = " << epsilon_D_c_cgs*(pow(fm_cgs, 3.0)*pow(c_cgs, 2.0))/MeV_cgs << endl;


    // RKF45 tolerances
    RKF45_tolerance_B = P_B_c*RKF45_tolerance_factor; // baryonic RKF45 tolerance
    RKF45_tolerance_D = P_D_c*RKF45_tolerance_factor; // DM RKF45 tolerance


    // dimensionless EOS file data
    vector<double> log10epsilon_B_EOSfile;
    vector<double> log10P_B_EOSfile;
    vector<double> log10h_B_EOSfile;
    vector<double> log10n_B_EOSfile;


    // convert log10epsilon and log10P from cgs to dimensionless
    for(int i = 0; i < log10epsilon_B_EOSfile_cgs.size(); i++)
    {
        log10epsilon_B_EOSfile.push_back(log10epsilon_B_EOSfile_cgs[i] + log10(kappa_cgs*G_cgs/pow(c_cgs, 2.0)));
        log10P_B_EOSfile.push_back(log10P_B_EOSfile_cgs[i] + log10(kappa_cgs*G_cgs/pow(c_cgs, 4.0)));
        log10h_B_EOSfile.push_back(log10h_B_EOSfile_cgs[i] - log10(pow(c_cgs, 2.0)));
        log10n_B_EOSfile.push_back(log10n_B_EOSfile_cgs[i] + log10(pow(kappa_cgs, 3.0/2.0)));
    }
    
//    tk::spline log10epsilon_B_of_log10h_B(log10h_B_EOSfile, log10epsilon_B_EOSfile, tk::spline::cspline, true);


//    // find central particle densities
//    double n_B_c = n_B(log10h_B_EOSfile, log10n_B_EOSfile, h_B_c, &n_nearest);
//    double n_D_c = n_D(P_D_c, m_D, y);


    // TOV solution vectors
    vector<double> r2;
    vector<double> h_B;
    vector<double> M_B;
    vector<double> h_D;
    vector<double> M_D;


    // TOV central conditions
    r2.push_back(r2_c);
    h_B.push_back(h_B_c);
    M_B.push_back(M_B_c);
    h_D.push_back(h_D_c);
    M_D.push_back(M_D_c);

    // solve TOV equations
    RK4_TOV(log10epsilon_B_EOSfile, log10P_B_EOSfile, log10h_B_EOSfile, r2, h_B, M_B, h_D, M_D, h_B_c, h_D_c, m_D, y);



    // find R_B and R_D
    int R_B_found = 0;
    int R_B_index;
    int R_D_found = 0;
    int R_D_index;
    for(int i = 0; i < r2.size(); i++)
    {
        if(h_B[i] < h_B_surf && R_B_found == 0)
        {
            R_B_index = i;
            R_B_found += 1;
        }
        if(h_D[i] <= h_D_surf_zero && R_D_found == 0)
        {
            R_D_index = i;
            R_D_found += 1;
        }
    }
    if(R_B_found == 0)
    {
        R_B_index = r2.size() - 1;
    }
    if(R_D_found == 0)
    {
        R_D_index = r2.size() - 1;
    }


    double R_B = sqrt(r2[R_B_index]); // R_B
    double R_D = sqrt(r2[R_D_index]); // R_D

    double M_B_of_R_B = M_B[R_B_index];          // M_B(R_B)
    double M_D_of_R_B = M_D[R_B_index];          // M_D(R_B)
    double M_T_of_R_B = M_B_of_R_B + M_D_of_R_B; // M_T(R_B)
    double M_T_of_R_B_over_R_B = M_T_of_R_B/R_B; // M_T(R_B)/R_B

    double M_B_of_R_D = M_B[R_D_index];          // M_B(R_D)
    double M_D_of_R_D = M_D[R_D_index];          // M_D(R_D)
    double M_T_of_R_D = M_B_of_R_D + M_D_of_R_D; // M_T(R_D)

    double M_T = M_B_of_R_B + M_D_of_R_D; // M_T
    double M_halo = M_T - M_T_of_R_B;     // M_halo
    double f_D = M_D_of_R_D/M_T;          // f_D

    // M_halo/R_D
    double M_halo_over_R_D;
    double M_halo_over_R_D_minus_R_B;
    double M_halo_over_M_D_of_R_D;
    if(R_D > R_B)
    {
        M_halo_over_R_D = M_halo/R_D;
        M_halo_over_R_D_minus_R_B = M_halo/(R_D - R_B);
        M_halo_over_M_D_of_R_D = M_halo/M_D_of_R_D;
    } else
    {
        M_halo_over_R_D = M_halo_over_R_D_zero;
        M_halo_over_R_D_minus_R_B = M_halo_over_R_D_minus_R_B_zero;
        M_halo_over_M_D_of_R_D = M_halo_over_M_D_of_R_D_zero;
    }



//    // metric vectors
//    vector<double> Phi;
//    vector<double> g;
//    vector<double> f;
//
//    // surface metric values
//    Phi.push_back(log(1.0 - 2.0*M_T/max(R_B, R_D))/2.0);
//    g.push_back(exp(2.0*Phi.back()));
//    f.push_back(1.0/g.back());

    // calculate metric throughout the star
    vector<double> r;
    for(int i = 0; i < r2.size(); i++)
    {
        r.push_back(sqrt(r2[i]));
    }
//    metric_integrate_trapezoid(r, P_B, M_B, P_D, M_D, Phi, g, f, epsilon_B_c, epsilon_D_c, R_B_index, R_D_index);
//
//    double g_R_B = g[R_B_index];                 // g(R_B)
//    double g_Sch_R_B = 1.0 - 2.0*M_T_of_R_B/R_B; // g_Sch(R_B)
//
//    vector<double> Phi_Sch;
//    vector<double> g_Sch;
//    vector<double> f_Sch;
//
//    if(R_B >= R_D)
//    {
//        for(int i = 0; i < g.size(); i++)
//        {
//            g_Sch.push_back(g[i]);
//            f_Sch.push_back(f[i]);
//        }
//    } else
//    {
//        g_Sch.push_back(g_Sch_R_B);
//        f_Sch.push_back(1.0/g_Sch_R_B);
//        Phi_Sch.push_back(log(g_Sch_R_B)/2.0);
//
//        metric_integrate_trapezoid(r, P_B, M_B, P_D, M_D, Phi_Sch, g_Sch, f_Sch, epsilon_B_c, epsilon_D_c, R_B_index, R_D_index_zero);
//
//        for(int i = (int)g_Sch.size(); i < g.size(); i++)
//        {
//            g_Sch.push_back(1.0 - 2.0*M_T_of_R_B/r[i]);
//            f_Sch.push_back(1.0/g_Sch.back());
//        }
//    }
//
//
//    double gravity_R_B = (1.0/sqrt(f[R_B_index]))*dPhidr(r[R_B_index], P_B[R_B_index], M_B[R_B_index], P_D[R_B_index], M_D[R_B_index], epsilon_B_c, P_B_c, epsilon_D_c, P_D_c); // gravity(R_B)
//    double gravity_Sch_R_B = M_T_of_R_B/(pow(R_B, 2.0)*sqrt(1.0 - 2.0*M_T_of_R_B/R_B));                                                                                        // gravity_Sch(R_B)
//
//
//    // particle number vectors
//    vector<double> N_B;
//    vector<double> N_D;
//
//    // particle number central conditions
//    N_B.push_back(N_B_c);
//    N_D.push_back(N_D_c);
//
//    // calculate particle numbers throughout the star
//    N_B_integrate_trapezoid(log10P_B_EOSfile, log10n_B_EOSfile, r, P_B, M_B, N_B, M_D, epsilon_B_c, n_B_c, epsilon_D_c, R_B_index);
//    N_D_integrate_trapezoid(r, M_B, P_D, M_D, N_D, m_D, y, epsilon_B_c, epsilon_D_c, n_D_c, R_D_index);
//
//    if(N_B.size() != N_D.size())
//    {
//        if(N_B.size() < N_D.size())
//        {
//            while(N_B.size() != N_D.size())
//            {
//                N_B.push_back(N_B.back());
//            }
//        } else
//        {
//            while(N_B.size() != N_D.size())
//            {
//                N_D.push_back(N_D.back());
//            }
//        }
//    }
//    double N_B_T = N_B[R_B_index]; // N_B(R_B)
//    double N_D_T = N_D[R_D_index]; // N_D(R_D)
//
//
//
//    vector<double> DANS_data = {R_B, R_D, M_B_of_R_B, M_T_of_R_B, M_D_of_R_D, M_T_of_R_D, M_T, M_halo, f_D, M_halo_over_R_D, M_halo_over_R_D_minus_R_B, M_halo_over_M_D_of_R_D, g_R_B, g_Sch_R_B, gravity_R_B, gravity_Sch_R_B, N_B_T, N_D_T};



    cout << "*******************************************" << endl;
    cout << "m_D = " << m_D_MeVfm << " MeV" << endl;
    cout << "y = " << y << endl;
    cout << "R_B = " << R_B*sqrt(kappa_cgs)/km_cgs << " km" << endl;
    cout << "R_D = " << R_D*sqrt(kappa_cgs)/km_cgs << " km" << endl;
    cout << "M_B(R_B) = " << M_B_of_R_B*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "M_T(R_B) = " << M_T_of_R_B*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "M_D(R_D) = " << M_D_of_R_D*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "M_T(R_D) = " << M_T_of_R_D*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "M_T = " << M_T*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "M_halo = " << M_halo*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << " solar masses" << endl;
    cout << "f_D = " << f_D << endl;
    cout << "M_halo/R_D = " << M_halo_over_R_D << endl;
    cout << "M_halo/(R_D - R_B) = " << M_halo_over_R_D_minus_R_B << endl;
    cout << "M_halo/M_D(R_D) = " << M_halo_over_M_D_of_R_D << endl;
    cout << "M_T(R_B)/R_B = " << M_T_of_R_B_over_R_B << endl;
//    cout << "g(R_B) = " << g_R_B << endl;
//    cout << "g_Sch(R_B) = " << g_Sch_R_B << endl;
//    cout << "Delta_rel g(R_B) = " << abs(g_R_B - g_Sch_R_B)/g_Sch_R_B << endl;
//    cout << "gravity(R_B) = " << gravity_R_B*pow(c_cgs, 2.0)/sqrt(kappa_cgs) << " cm/s^2" << endl;
//    cout << "gravity_Sch(R_B) = " << gravity_Sch_R_B*pow(c_cgs, 2.0)/sqrt(kappa_cgs) << " cm/s^2" << endl;
//    cout << "Delta_rel gravity(R_B) = " << abs(gravity_R_B - gravity_Sch_R_B)/gravity_Sch_R_B << endl;
//    cout << "N_B_T = " << N_B_T << endl;
//    cout << "N_D_T = " << N_D_T << endl;
    cout << "*******************************************" << endl;

    
    
    // save DANS data in file
    remove(DANS_data_filename.c_str());
    ofstream DANS_data_file;
    DANS_data_file.open(DANS_data_filename);
//    for(int i = 0; i < r.size(); i++)
//    for(int i = 0; i <= max(R_B_index, R_D_index); i++)
//    {
//        DANS_data_file << setprecision(12) << r[i]*sqrt(kappa_cgs)/km_cgs << "\t" << epsilon_B(log10epsilon_B_EOSfile, log10P_B_EOSfile, P_B[i])*pow(c_cgs, 2.0)/(kappa_cgs*G_cgs)*pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)/MeV_cgs << "\t" << P_B[i]*pow(c_cgs, 4.0)/(kappa_cgs*G_cgs)*pow(fm_cgs, 3.0)/MeV_cgs << "\t" << M_B[i]*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << "\t" << N_B[i] << "\t" << epsilon_D(P_D[i], m_D, y, "dimensionless")*pow(c_cgs, 2.0)/(kappa_cgs*G_cgs)*pow(fm_cgs, 3.0)*pow(c_cgs, 2.0)/MeV_cgs << "\t" << P_D[i]*pow(c_cgs, 4.0)/(kappa_cgs*G_cgs)*pow(fm_cgs, 3.0)/MeV_cgs << "\t" << M_D[i]*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << "\t" << N_D[i] << "\t" << g[i] << "\t" << g_Sch[i] << "\t" << f[i] << "\t" << f_Sch[i] << endl;
//    }
    for(int i = 0; i <= max(R_B_index, R_D_index); i++)
    {
        DANS_data_file << setprecision(12)
        << r[i]*sqrt(kappa_cgs)/km_cgs << "\t"
        << epsilon_B(log10epsilon_B_EOSfile, log10h_B_EOSfile, h_B[i], &n_nearest)*pow(c_cgs, 2.0)/(kappa_cgs*G_cgs) << "\t"
        << P_B(log10P_B_EOSfile, log10h_B_EOSfile, h_B[i], &n_nearest)*pow(c_cgs, 4.0)/(kappa_cgs*G_cgs) << "\t"
        << h_B[i]*pow(c_cgs, 2.0) << "\t"
        << M_B[i]*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << "\t"
        << 0.0 << "\t"// << N_B[i] << "\t"
        << epsilon_D(h_D[i], m_D, y, "dimensionless")*pow(c_cgs, 2.0)/(kappa_cgs*G_cgs) << "\t"
        << P_D(h_D[i], m_D, y, "dimensionless")*pow(c_cgs, 4.0)/(kappa_cgs*G_cgs) << "\t"
        << h_D[i]*pow(c_cgs, 2.0) << "\t"
        << M_D[i]*sqrt(kappa_cgs)*pow(c_cgs, 2.0)/(G_cgs*solarmass_cgs) << "\t"
        << 0.0 << "\t"// << N_D[i] << "\t"
        << 0.0 << "\t"// << g[i] << "\t"
        << 0.0 << "\t"// << g_Sch[i] << "\t"
        << 0.0 << "\t"// << f[i] << "\t"
        << 0.0 << endl;// << f_Sch[i] << endl;
    }
    DANS_data_file.close();
    
//    ofstream enthalpy_data_file;
//    enthalpy_data_file.open(enthalpy_data_filename, std::ios_base::app);
//    enthalpy_data_file << setprecision(6) << h_D_c/h_B_c << "\t" << R_D/R_B << endl;
//    enthalpy_data_file.close();



//    return DANS_data;
}




/***************************************************************************/
/* Routine that locates nearest grid point for a given value.              */
/* Adapted from Numerical Recipes.                                         */
/***************************************************************************/
void hunt(const vector<double>& xx, int n, double x, int *jlo)
{
    int jm,jhi,inc,ascnd;

    ascnd=(xx[n - 1] >= xx[0]);
    if (*jlo < 0 || *jlo > n - 1) {
        *jlo=0;
        jhi=n-1;
    } else {
        inc=1;
        if (x >= xx[*jlo] == ascnd) {
            if (*jlo == n - 1) return;
            jhi=(*jlo) + 1;
            while (x >= xx[jhi] == ascnd) {
                *jlo=jhi;
                inc += inc;
                jhi=(*jlo)+inc;
                if (jhi > n - 1) {
                    jhi=n - 1;
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
            while (x < xx[*jlo] == ascnd) {
                jhi=(*jlo);
                inc += inc;
                *jlo=jhi-inc;
                if (*jlo < 0) {
                    *jlo=0;
                    break;
                }
            }
        }
    }
    while (jhi-(*jlo) > 1) {
        jm=(jhi+(*jlo)) >> 1;
        if (x > xx[jm] == ascnd)
            *jlo=jm;
        else
            jhi=jm;
    }
}

/*C*/
/*************************************************************************/
/* Driver for the interpolation routine. First we find the tab. point    */
/* nearest to xb, then we interpolate using four points around xb.       */
/*************************************************************************/
double interp(const vector<double>& xp,
              const vector<double>& yp,
              int    np ,
              double xb,
              int    *n_nearest_pt)
{

 int k,        /* index of 1st point */
     m=4;      /* degree of interpolation */
 
 double y;     /* intermediate value */

 hunt(xp,np,xb,n_nearest_pt);
//    cout << "Hunt complete" << endl;
 k=IMIN(IMAX((*n_nearest_pt)-(m-1)/2,0),np-m);

 if( xb==xp[k] ||  xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3])
 {
     xb += DBL_EPSILON;
 }

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
