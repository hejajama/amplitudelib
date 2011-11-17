/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

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
#include <gsl/gsl_min.h>


#include <algorithm>


extern "C"
{
    #include "../fourier/fourier.h"
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 *
 * By default der=0, der=1 is 1st derivative w.r.t r, 2 2nd
 * By default bspline=false, it it is true use bspline interpolation
 *  (useful for noisy data), affects only if there is no previously allocated
 *  interpolator
 */
REAL AmplitudeLib::N(REAL r, REAL y, int der, bool bspline)
{    
    if (der>2 or der<0)
    {
        cerr << "Derivative degree " << der << " is not 0, 1 or 2! " << LINEINFO
         << endl;
         return 0;
    }
    if (r < MinR() or r > MaxR() )
    {
        if (out_of_range_errors)
            cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
                << " asked r=" << r << " " << LINEINFO
                << endl;
        if (r<MinR()) r=MinR()*1.000001; else if (r>MaxR()) r=MaxR()*0.999999;
    }

    if (y<0 or y>yvals[yvals.size()-1] )
    {
        if (out_of_range_errors)
            cerr << "y must be between limits [" << 0 << ", "
                << yvals[yvals.size()-1] << "], asked y=" << y << " "
                << LINEINFO << endl;
        if (y < 0) y=0; else if (y>yvals[yvals.size()-1]) y=yvals[yvals.size()-1];
    }
    //REAL lnr = std::log(r);

    /// Use already initialized interpolator
    if (std::abs(y - interpolator_y) < 0.01 and !bspline )
    {
        REAL result=0;
        if (der==0) 
            result = interpolator->Evaluate(r);
        if (der==1)
        {
            result = interpolator->Derivative(r);
        }
        if (der==2)
        {
            result = interpolator->Derivative2(r);
        }
        if (!kspace and result>1 and der==0) return 1;  // in x space limit N(r)<=1
        if (result<0 and der==0) return 0;              // limit N >= 0
        return result;

    }

    /// Initialize new interpolator and use it
    int yind = FindIndex(y, yvals);
    int rind = FindIndex(r, rvals);

    int interpolation_points = INTERPOLATION_POINTS;
    if (bspline) interpolation_points += 10;

    int interpolation_start, interpolation_end;
    if (rind - interpolation_points/2 < 0)
    {
		interpolation_start=0;
		interpolation_end=interpolation_points;
	}
	else if (rind + interpolation_points/2 > RPoints()-1 )
	{
		interpolation_end = RPoints()-1;
		interpolation_start = RPoints()-interpolation_points-3;
	}
	else
	{
		interpolation_start = rind - interpolation_points/2;
		interpolation_end = rind + interpolation_points/2;
	}

    int interpo_points = interpolation_end - interpolation_start+1;
    REAL *tmparray = new REAL[interpo_points];
    REAL *tmpxarray = new REAL[interpo_points];
    for (int i=interpolation_start; i<= interpolation_end; i++)
    {
		tmpxarray[i-interpolation_start]=rvals[i];

        tmparray[i-interpolation_start] = n[yind][i];

        // Interpolate in y if possible
		if (yind < static_cast<int>(yvals.size()-1) )
        {
            tmparray[i-interpolation_start]=n[yind][i] 
                + (y - yvals[yind]) * (n[yind+1][i] - n[yind][i])
                / (yvals[yind+1]-yvals[yind]);
		}
    }
    
    Interpolator interp(tmpxarray, tmparray, interpo_points);
    if (bspline) interp.SetMethod(INTERPOLATE_BSPLINE);
    interp.Initialize();
    REAL result=0;
    if (der==0)
        result = interp.Evaluate(r);
    if (der==1)
    {
        result = interp.Derivative(r);
        //result /= r;    // dN / d ln r = r dN/dr
    }
    if (der==2)
    {
        result = interp.Derivative2(r);
        //result /= SQR(r);   // d^2N / d lnr^2 = r^2 d^2 N / dr^2
    }
    
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

REAL AmplitudeLib::S(REAL r, REAL y, int der)
{
    double s = 1.0 - N(r,y,der);
    if (s<=1e-14) return 0;
    return s;
}

/*
 * FT the amplitude to the k-space
 * N(k) = \int d^2 r/(2\pi r^2) exp(ik.r)N(r)
 *  = \int dr/r BesselJ[0,k*r] * N(r)
 * 
 * Note: for performance reasons it is probably a good idea to
 * call AmplitudeLib::InitializeInterpolation(y) before this
 */
struct N_k_helper
{
    REAL y; REAL kt;
    bool adjoint;
    AmplitudeLib* N;
};
REAL N_k_helperf(REAL r, void* p);
REAL AmplitudeLib::N_k(REAL kt, REAL y)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    SetOutOfRangeErrors(false);
    
    N_k_helper par;
    par.y=y; par.N=this; par.kt=kt;
    REAL result = fourier_j0(kt,N_k_helperf,&par);
    return result;
}

REAL N_k_helperf(REAL r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    if (r < par->N->MinR()) return 0;
    else if (r > par->N->MaxR()) return 1.0/r;
    return 1.0/r*par->N->N(r, par->y);
}

/*
 * FT amplitude from k to x space
 * N(r) = r^2 \int d^2 k/(2\pi) exp(-ik.r) N(k)
 *      = r^2 \int dk k BesselJ[0,k*r] * N(k)
 */
struct N_k_to_x_helper
{
    REAL y; REAL x;
    AmplitudeLib* N;
};
REAL N_k_to_x_helperf(REAL k, void* p)
{
    N_k_to_x_helper* par = (N_k_to_x_helper*)p;
    return par->N->N(k, par->y)*k;
}
REAL AmplitudeLib::N_k_to_x(REAL x, REAL y)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    SetOutOfRangeErrors(false);
    
    N_k_to_x_helper par;
    par.y=y; par.N=this; par.x=x;
    REAL result = fourier_j0(x,N_k_to_x_helperf,&par);
    return x*x*result;
}


/*
 * FT S=1-N to the k-space
 * S(k) = \int d^2 r exp(ik.r) (1-N(r))
 *  = (2\pi) \int dr r BesselJ[0,k*r] * (1-N(r))
 * 
 * Note: for performance reasons it is probably a good idea to
 * call AmplitudeLib::InitializeInterpolation(y) before this
 *
 * if adjoint==true, then use adjoint representation for N
 * Default value for adjoint is false
 */
REAL S_k_helperf(REAL r, void* p);
REAL AmplitudeLib::S_k(REAL kt, REAL y, bool adjoint)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    
    N_k_helper par;
    par.y=y; par.N=this; par.kt=kt; par.adjoint=adjoint;

    REAL result=0;

    if (kt < 1e-5)  // k_T \approx 0 -> integrate just \int d^2 r S(r)
    {
        gsl_function fun; fun.function=S_k_helperf;
        fun.params=&par;
    

        REAL abserr; size_t eval;
        int status = gsl_integration_qng(&fun, MinR(), MaxR(),
            0, 0.01,  &result, &abserr, &eval);
        if (status)
        {
            cerr << "S_k integral failed at " << LINEINFO <<", y=" << y
            << " k_T=" << kt << ", adjoint: " << adjoint << endl;
        }
    }
    else
        result = fourier_j0(kt,S_k_helperf,&par);
    return result*2.0*M_PI; 
}

REAL S_k_helperf(REAL r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    REAL result;
    if (!par->adjoint)
    {
        if (r > par->N->MaxR()) return 0;
        result = r*(1.0-par->N->N(r, par->y));
    }
    else
    {
        result = r*(1.0-par->N->N_A(r, par->y));
    }

    return result;
}

/*
 * Virtual photon-proton cross sections
 * Separately for transversially and longitudinally polarized photons
 */
struct Inthelper_totxs
{
    AmplitudeLib* N;
    int pol;    // 0: L, 1: T
    REAL Qsqr,y;
    WaveFunction* wf;
};

REAL Inthelperf_totxs(REAL r, void* p)
{
    Inthelper_totxs* par = (Inthelper_totxs*)p;

    REAL result = r*par->N->N(r,par->y);
    if (par->pol==0)    // Longitudinal
        result *= par->wf->PsiSqr_L_intz(par->Qsqr, r);
    else if (par->pol==1)   // Transverse
        result *= par->wf->PsiSqr_T_intz(par->Qsqr, r);
    else
        cerr << "Invalid polarization " << par->pol << " at " << LINEINFO << endl;

    return result;
}

REAL AmplitudeLib::ProtonPhotonCrossSection(REAL Qsqr, REAL y, int pol)
{
    Inthelper_totxs par; par.N=this;
    par.pol=pol; par.Qsqr=Qsqr; par.y=y;
    VirtualPhoton wavef;
    par.wf=&wavef;

    gsl_function fun; fun.function=Inthelperf_totxs;
    fun.params=&par;
    

    REAL result,abserr; size_t eval;
    int status = gsl_integration_qng(&fun, MinR(), MaxR(),
    0, 0.01,  &result, &abserr, &eval);
    
    if(status){ std::cerr<< "r integral in ProtonPhotonCrossSection failed with code " 
        << status << " (Qsqr=" << Qsqr << ", y=" << y 
        << " relerr=" << abserr/result << ") at " << LINEINFO << std::endl;
    }

    return 2.0*M_PI*result; //2\pi from \theta integral

    

}


/*
 * d ln N / d ln r^2 = 1/N * d N / d ln r^2 = 1/N d N / dr^2  dr^2 / d ln r^2
 *  = 1/N d N / dr^2 r^2 = 1/N d N / dr  dr / dr^2  r^2 = r/(2N) * dN/dr
 */
REAL AmplitudeLib::LogLogDerivative(REAL r, REAL y)
{
    REAL dndr = N(r,y,1);
    return r/(2.0*N(r,y))*dndr;
}

/*
 * Calculate saturation scale, definition is
 * N(r_s, y) = Ns = (for example) 0.5
 *
 * In kspace: satscale = max of k^(2\gamma_c), \gamma_c=0.6275
 * Ref. hep-ph/0504080
 */
struct SatscaleSolverHelper
{
    REAL y;
    REAL Ns;
    REAL gammac;
    AmplitudeLib* N;
};
REAL SatscaleSolverHelperf(double r, void* p)
{
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    return par->N->N(r, par->y) - par->Ns;
}

REAL SatscaleSolverHelperf_k(REAL kt, void* p)
{
   
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    REAL val = -std::pow(kt, 2.0*par->gammac)
        *par->N->N(kt, par->y);
    return val;
    
}


REAL AmplitudeLib::SaturationScale(REAL y, REAL Ns)
{
    
    SatscaleSolverHelper helper;
    helper.y=y; helper.Ns=Ns; helper.N=this; helper.gammac=0.6275;
    const int MAX_ITER = 300;
    const REAL ROOTFINDACCURACY = 0.01;
    gsl_function f;
    f.params = &helper;
        
    f.function = &SatscaleSolverHelperf;

    const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        
    gsl_root_fsolver_set(s, &f, MinR()*1.0001, MaxR()*0.999);
    int iter=0; int status; double min,max;
    do
    {
        iter++;
        gsl_root_fsolver_iterate(s);
        min = gsl_root_fsolver_x_lower(s);
        max = gsl_root_fsolver_x_upper(s);
        status = gsl_root_test_interval(min, max, 0, ROOTFINDACCURACY);    
    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (iter>=MAX_ITER)
        cerr << "Solving failed at y=" << y << " " << LINEINFO << endl;


    REAL sat = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    if (!kspace)
        return sat;
    

    
    f.function=&SatscaleSolverHelperf_k;
  
    const gsl_min_fminimizer_type *Tm = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *sm = gsl_min_fminimizer_alloc(Tm);

    // Interpolation doesn't work well at the edge of the ktsqr range, so limit
    // the studied interval
    REAL interval_min = MinR()*1.0001;
    REAL interval_max = MaxR()*0.999;

    //gsl_error_handler_t* handler = gsl_set_error_handler_off();
    if (gsl_min_fminimizer_set(sm, &f, sat, interval_min, interval_max)
        == GSL_EINVAL)
    {
        // Endpoints do not enclose a minimum (on the other hand
        // helper(pos) > helper(min),helper(max), but we can anyway continue
        cerr << "Endpoints do not enclose a minumum! " << LINEINFO << endl;
    }
    //gsl_set_error_handler(handler);


    iter=0;
    do
    {
        iter++;
        gsl_min_fminimizer_iterate(sm);
        sat = gsl_min_fminimizer_x_minimum(sm);
        interval_min = gsl_min_fminimizer_x_lower (sm);
        interval_max = gsl_min_fminimizer_x_upper (sm);
        ///TODO: relerror/abserror from ktsqrval difference
        status = gsl_min_test_interval (interval_min, interval_max, 0.0, 0.01);
        
    } while (status == GSL_CONTINUE and iter < MAX_ITER);

    if (status == GSL_CONTINUE)
    {
        cerr << "Didn't find saturation scale when y=" << y  << ", iterated "
            << iter << "/" << MAX_ITER << " times at " << LINEINFO << endl;
    }

    gsl_min_fminimizer_free(sm);

    //cout << "start: " << std::exp(0.4*y) << " min: " << Ktsqrval(1) << " max: " << Ktsqrval(KtsqrPoints()-2)
    //    << " minimumsqr: " << pos << " iterations: " << iter << endl;

    return sat;
}

/*
 * Amplitude in adjoint representation
 * Coordinate space: N_A(r) = 2N(r)-N(r)^2
 */
REAL AmplitudeLib::N_A(REAL r, REAL y, int der)
{
    if (kspace)
    {
        cerr <<"N_A is not implemented for mome. space!" << endl;
        return 0;
    }
    if (der<0 or der>2)
    {
        cerr << "Derivative order " << der << " is not valid. " << LINEINFO << endl;
        return -1;
    }
    REAL result=0;

    if (der==0)
    {
        REAL n = N(r,y);
        result = 2.0*n - n*n;
    }
    else if (der==1)
    {
        REAL n = N(r,y);
        REAL dn = N(r,y,1);
        result= 2.0*dn - 2.0*n*dn;
    }
    else if (der==2)
    {
        REAL n = N(r,y);
        REAL dn = N(r,y,1);
        REAL d2n = N(r,y,2);
        result = 2.0*d2n - 2.0*dn*dn - 2.0*n*d2n;
    }
    

    return result;
}


/*
 * Load data from a given file
 * Format is specified in file bk/README
 * If kspace is true, the datafile is in kspace solved using BK code
 * ("x axis" in kt^2 in datafile, don't limit N<=1)
 */
AmplitudeLib::AmplitudeLib(std::string datafile, bool kspace_)
{
    out_of_range_errors=true;
    kspace=kspace_;
    DataFile data(datafile);
    data.GetData(n, yvals);
    minr = data.MinR();
    rmultiplier = data.RMultiplier();
    rpoints = data.RPoints();
    x0 = data.X0();

    tmprarray = new REAL[data.RPoints()];
    tmpnarray = new REAL[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        REAL tmpr = std::log(minr*std::pow(rmultiplier, i));
        if (kspace)
            tmpr = 0.5*tmpr;    // take square root if kspace, so k^2 -> k
        lnrvals.push_back(tmpr);
        rvals.push_back(std::exp(lnrvals[i]));
        tmprarray[i] = std::exp(lnrvals[i]);
    }

    if (kspace)
    {
        minr=std::sqrt(minr);
        rmultiplier = std::sqrt(rmultiplier);
    }
    
    cout << "# Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << " x0 " << X0() << endl;

    interpolator_y=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)
}

/*
 * Initializes interpolation method with all read data points at given y
 * If bspline is true, then use bspline interpolation, default is false
 */
void AmplitudeLib::InitializeInterpolation(REAL y, bool bspline)
{
    if (std::abs(interpolator_y - y) < 0.01) return;    // Already done it
    if (interpolator_y>=0)
    {
        interpolator_y=-1;
        delete interpolator;
    }
    for (int i=0; i<rpoints; i++)
    {
        REAL tmpr = tmprarray[i];
        if (i==0) tmpr*=1.001; if (i==rpoints-1) tmpr*=0.999;
        tmpnarray[i] = N(tmpr, y);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);
    if (bspline)
        interpolator->SetMethod(INTERPOLATE_BSPLINE);
    interpolator->Initialize();
    interpolator_y = y;
}

/*
 * Initializes interpolator and returns it
 */
Interpolator* AmplitudeLib::MakeInterpolator(REAL y)
{
    REAL* interp_narray = new REAL[rpoints];
    for (int i=0; i<rpoints; i++)
    {
        tmpnarray[i] = N(tmprarray[i], y);
    }
    Interpolator* inter = new Interpolator(tmprarray, tmpnarray, rpoints);
    inter->Initialize();
    delete[] interp_narray;
    return inter;
        
}

AmplitudeLib::~AmplitudeLib()
{
    if (interpolator_y>=0)
    {
        delete interpolator;
    }
    delete[] tmpnarray;
    delete[] tmprarray;
}

int AmplitudeLib::RPoints()
{
    return rpoints;
}
REAL AmplitudeLib::MinR()
{
    return minr;
}

REAL AmplitudeLib::MaxR()
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int AmplitudeLib::YVals()
{
    return yvals.size();
}

REAL AmplitudeLib::MaxY()
{
    return yvals[yvals.size()-1];
}

void AmplitudeLib::SetOutOfRangeErrors(bool er)
{
    out_of_range_errors=er;
}

REAL AmplitudeLib::X0()
{
    return x0;
}

