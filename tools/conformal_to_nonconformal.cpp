/*
 * Calculate non-conformal dipole from conformal dipole
 * (up to as^2 corrections)
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
 */

#include <iostream>
#include "../amplitudelib/amplitudelib.hpp"
#include <string>
#include <cmath>
#include <gsl/gsl_integration.h>
#include "tools.hpp"
#include "../amplitudelib/qcd.hpp"
#include "tools.hpp"

using namespace Amplitude;
using namespace std;

const double a = 10;

struct Inthelper
{
    AmplitudeLib* N;    // conformal dipole
    double r;   // = x - y = parent dipole
    double xbj;
    double z;   // = y - z = daughter dipole
};

double Inthelperf_lo_z(double v, void* p);
double Inthelperf_lo_theta(double theta, void* p);

const int RINTPOINTS=20;
const int THETAINTPOINTS=10;
const double INTACCURACY = 0.01;

QCD qcd;

double NonConfIntegral(double r, double x, AmplitudeLib* N)
{
    gsl_function fun;
    Inthelper helper;
    helper.r=r;
    helper.xbj=x;
    helper.N=N;
    
    fun.params = &helper;
    fun.function = Inthelperf_lo_z;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(RINTPOINTS);

    double minlnr = std::log( 0.5*N->MinR() );
    double maxlnr = std::log( 2.0*N->MaxR() );

    int status; double  result, abserr;
    status=gsl_integration_qag(&fun, minlnr,
            maxlnr, 0, INTACCURACY, RINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //cerr << "RInt failed, r=" << r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //    << std::abs(abserr/result) << endl;
    }

    
    double s = N->S(r, x) + result;
    return 1.0-s;
}

double Inthelperf_lo_z(double z, void* p)
{
    Inthelper* helper = reinterpret_cast<Inthelper*>(p);
    helper->z=std::exp(z);
    //cout << " r " << helper->r << " v " << helper->v << endl;
    gsl_function fun;
    fun.function=Inthelperf_lo_theta;
    fun.params = helper;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0,
            2.0*M_PI, 0, INTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS21, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        //#pragma omp critical
        //cerr << "RInt failed, v=" << v <<", r=" << helper->r <<": at " << LINEINFO << ", result " << result << ", relerr "
        //<< std::abs(abserr/result) << endl;
    }

    result *= std::exp(2.0*z);  // Jacobian v^2 dv
    //result *= 2.0;  // As integration limits are just [0,pi]
    return result;
}

double Inthelperf_lo_theta(double theta, void* p)
{
        //cout << " theta " << theta << endl;
    Inthelper* par = reinterpret_cast<Inthelper*>(p);

    double result=0;
    double r = par->r;
    double z = par->z;
    double xbj = par->xbj;

    double r_m_z = sqrt( r*r + z*z - 2.0*r*z*cos(theta));

    double musqr = 4.0/(r*r);
    result = qcd.Alphas(musqr)*qcd.Nc()/(4.0*M_PI*M_PI) * SQR(r/(z*r_m_z)) * log( SQR( a*r/(z*r_m_z) ) )
            * ( par->N->S(z, xbj) * par->N->S(r_m_z, xbj) - par->N->S(r, xbj) );
    return result;

     
}

int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    string fname = argv[1];
    AmplitudeLib N(fname);



    N.SetOutOfRangeErrors(false);
    
    double y = StrToReal(argv[2]);
    double xbj = N.X0()*exp(-y);

    for (double r = N.MinR(); r<N.MaxR(); r*=1.1)
    {
        cout << r << " " << NonConfIntegral(r, xbj, &N) << endl;
    }
    return 0;

}
