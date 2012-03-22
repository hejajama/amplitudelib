#include "amplitudelib.hpp"
#include "virtual_photon.hpp"
#include "datafile.hpp"
#include "../tools/tools.hpp"
#include "../tools/config.hpp"
#include <string>
#include <cmath>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_sf_bessel.h>
extern "C"
{
    #include "../fourier/fourier.h"
}


/*
 * Unintegrated gluon density
 * \psi(k, y) = C_F/( (2\pi)^3 \alphas(k) ) * \int d^2 r exp(-ik.r) \nabla^2 N_G,
 * \nabla^2 N_G = (1/r \der_r + \der^2_r)(2N-N^2)
 *  = 2/r \der N + 2 \der^2 N - 2N/r \der N - 2(\der r N)^2 - 2N \der^2 N
 *
 * For performance reasons it is probably a good idea to call
 * InitializeInterpoaltion before this
 */

struct UGDHelper
{
    double y;
    AmplitudeLib* N;
    Interpolator* interp;
};
double UGDHelperf(double x, void* p);

/*
 * Calculate unintegrated gluon distribution
 * def. as C_f/\alpha_s(k) 1/(2\pi)^3 \int dr r \nabla_r^2 (2N-N^2)
 *
 * If interp != null, we know that it is already initialized with given y
 * and thus we can use it. Default value of interp is NULL
 */
double AmplitudeLib::UGD(double k, double y, Interpolator* interp)
{
    //if (k < UGD_IR_CUTOFF) return 0;
    
    set_fpu_state();
    init_workspace_fourier(700);   // number of bessel zeroes, max 2000

    UGDHelper par;
    par.y=y; par.N=this; par.interp=interp;
    double result = fourier_j0(k,UGDHelperf,&par);

    double Cf = (SQR(Nc)-1.0)/(2.0*Nc);
    result *= Cf/Alpha_s(SQR(k));  //Todo: scaling?
    result /= SQR(2.0*M_PI);

    return result;

}

double UGDHelperf(double r, void* p)
{
    UGDHelper* par = (UGDHelper*) p;
    double result=0;
    double dern=0, der2n=0, n=0;
    if (r >= par->N->MaxR()) // N==1 for all r>MaxR()
    {
        n=1; dern=0; der2n=0;
    }
    else if (r <= par->N->MinR())  
        return 0;
    else
    {
        if (par->interp == NULL)
        {
            dern = par->N->N(r, par->y, 1);
            der2n = par->N->N(r, par->y, 2);
            n = par->N->N(r, par->y);
        }
        else
        {
            dern = par->interp->Derivative(r);
            der2n = par->interp->Derivative2(r);
            n = par->interp->Evaluate(r);
        }
    }

    result = 2.0/r*dern + 2.0*der2n - 2.0*n/r*dern - 2.0*SQR(dern) - 2.0*n*der2n;
    return r*result;
}

/*
 * k_T factorization, calculate d\sigma^{A+B -> g} / (dy d^2 p_T)
 * Ref. 1011.5161 eq. (8)
 * Normalization is arbitrary
 */
struct Inthelper_dsigmadyd2pt
{
    AmplitudeLib* N;
    Interpolator* amplitude_y1;
    Interpolator* amplitude_y2;
    double pt;
    double y1,y2;
    double kt;
};
double Inthelperf_dN_gluon_dyd2pt_kint(double kt, void* p);
double Inthelperf_dN_gluon_dyd2pt_thetaint(double kt, void* p);
double AmplitudeLib::dN_gluon_dyd2pt(double pt, double y, double sqrts)
{
    const int KTINTPOINTS = 100;
    const double KTINTACCURACY = 0.05;
    
    Inthelper_dsigmadyd2pt helper;
    helper.N=this; helper.pt=pt;
    helper.amplitude_y1=NULL; helper.amplitude_y2=NULL;
    double y1 = pt/sqrts*std::exp(y); double y2 = pt/sqrts*std::exp(-y);
    helper.y1=y1; helper.y2=y2;
    if (y1<0 or y2<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", y1=" << y1 <<", y2="
            << y2 << endl;
        return 0;
    }

    helper.amplitude_y1 = MakeInterpolator(y1);
    helper.amplitude_y2 = MakeInterpolator(y2);

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_dN_gluon_dyd2pt_kint;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(KTINTPOINTS);

    double minkt = 0;
    double maxkt = pt;

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, minkt, maxkt,
            0, KTINTACCURACY, KTINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "k_T integration failed at " << LINEINFO <<", pt=" << pt
            << ", y1=" << y1 <<", y2=" << y2 <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }

    delete helper.amplitude_y1;
    delete helper.amplitude_y2;
    
    return result/(4.0*SQR(pt));    // 4.0 is in the integration measure d^2k/4
}

double Inthelperf_dN_gluon_dyd2pt_kint(double kt, void* p)
{
    const int THETAINTPOINTS = 40;
    const double THETAINTACCURACY = 0.05;
    
    Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*) p;
    par->kt=kt;
    gsl_function fun; fun.params=par;
    fun.function=Inthelperf_dN_gluon_dyd2pt_thetaint;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(THETAINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, 0, M_PI,
            0, THETAINTACCURACY, THETAINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "thetaintegration failed at " << LINEINFO <<", pt=" << par->pt
            << ", kt=" << kt <<", y1=" << par->y1 <<", y2=" << par->y2
            <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }

    return 2.0*result;  // 2.0 from int. limits [0,2\pi] -> [0,\pi]
    
}

double Inthelperf_dN_gluon_dyd2pt_thetaint(double theta, void* p)
{
    Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*) p;
    double costheta = std::cos(theta);
    double ktpluspt = SQR(par->kt) + SQR(par->pt) + 2.0*par->kt*par->pt*costheta;
    ktpluspt = std::sqrt(ktpluspt);
    double ktminuspt = SQR(par->kt) + SQR(par->pt) - 2.0*par->kt*par->pt*costheta;
    ktminuspt = std::sqrt(ktminuspt);
    
    double Q = std::max(ktpluspt/2.0, ktminuspt/2.0);
    double ugd1 = par->N->UGD(ktpluspt/2.0, par->y1, par->amplitude_y1);
    double ugd2 = par->N->UGD(ktminuspt/2.0, par->y2, par->amplitude_y2);

    if (isnan(ugd1) or isnan(ugd2) or Q<1e-10)
    {
        cerr << "???\n";
    }

    return Alpha_s(SQR(Q))*ugd1*ugd2;
}

/*
 * Rapidity distribution
 * d\sigma/dy = \int d^2 pt dSigmadydp2t
 */
double Inthelperf_dsigmady(double pt, void* p);
double AmplitudeLib::dSigmady(double y, double sqrts)
{
    /*double PTINTPOINTS = 400;
    double PTINTACCURACY = 0.05;
    
    double minpt = 0;
    double maxpt = 12;    // As in ref. 1011.5161

    Inthelper_dsigmadyd2pt helper;
    helper.N=this; helper.sqrts=sqrts; helper.y=y;

    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_dsigmadyd2pt;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(PTINTPOINTS);

    int status; double result, abserr;
    status=gsl_integration_qag(&fun, minpt, maxpt,
            0, PTINTACCURACY, PTINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
    {
        cerr << "p_T integration failed at " << LINEINFO 
            << ", y=" << y <<", result=" << result
            <<" relerr " << std::abs(abserr/result) << endl;
    }
*/
    return 0;
}

double Inthelperf_dsigmady(double pt, void* p)
{
    /*Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*)p;
    double x1 = pt/par->sqrts*std::exp(par->y);
    double x2 = pt/par->sqrts*std::exp(-par->y);

    return par->N->dN_gluon_dyd2pt(pt, x1, x2);
    */
    return 0;
}

/*
 * d\sigma/dy, monte carlo integration
 */

// x[0]=p, x[1]=k, x[2]=theta
double Inthelperf_dsigmadymc(double* x, size_t dim, void* p)
{
    /*Inthelper_dsigmadyd2pt* par = (Inthelper_dsigmadyd2pt*)p;
    if (x[1]>x[0])  return 0;
    

    double x1 = x[0]/par->sqrts*std::exp(par->y);
    double x2 = x[0]/par->sqrts*std::exp(-par->y);
    double costheta = std::cos(x[2]);
    double ktpluspt = SQR(x[1]) + SQR(x[0]) + 2.0*x[1]*x[0]*costheta;
    ktpluspt = std::sqrt(ktpluspt);
    double ktminuspt = SQR(x[1]) + SQR(x[0]) - 2.0*x[1]*x[0]*costheta;
    ktminuspt = std::sqrt(ktminuspt);
    
    double Q = std::max(ktpluspt/2.0, ktminuspt/2.0);
    double ugd1 = par->N->UGD(ktpluspt/2.0, x1);
    double ugd2 = par->N->UGD(ktminuspt/2.0, x2);

    return 1.0/x[0]*x[1]*Alpha_s(SQR(Q))*ugd1*ugd2;
    */
    return 0;
}
 
double AmplitudeLib::dSigmady_mc(double y, double sqrts)
{
    /*
    double maxpt = 12;    // As in ref. 1011.5161

    double lower[3] = {1e-8, 1e-8, 0};
    double upper[3] = {maxpt, maxpt, 2*M_PI};
    Inthelper_dsigmadyd2pt helper;
    helper.N=this; helper.y=y; helper.sqrts=sqrts;

    double res,err;
    
    const gsl_rng_type *T;
    gsl_rng *r;
     
    gsl_monte_function G = { &Inthelperf_dsigmadymc, 3, &helper };
     
    size_t calls = 500;
     
    gsl_rng_env_setup ();
     
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc (3);
    gsl_monte_miser_integrate (&G, lower, upper, 3, calls, r, s, 
                               &res, &err);
    gsl_monte_miser_free (s);

    return res;
    */
    return 0;

}
