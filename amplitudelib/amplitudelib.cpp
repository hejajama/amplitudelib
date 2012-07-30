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
double AmplitudeLib::N(double r, double y, int der, bool bspline)
{    
    if (bspline)
        cerr << "AmplitudeLib::N() with bspline is not supported (may not work) " << LINEINFO << endl;
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
        //if (out_of_range_errors)
            cerr << "y must be between limits [" << 0 << ", "
                << yvals[yvals.size()-1] << "], asked y=" << y << " "
                << LINEINFO << endl;
        if (y < 0) y=0; else if (y>yvals[yvals.size()-1]) y=yvals[yvals.size()-1];
    }
    //double lnr = std::log(r);

    
    /// Use already initialized interpolator
    if (std::abs(y - interpolator_y) < 0.01 and !bspline )
    {
        double result=0;
        if (der==0)
        { 
            if (r >= maxr_interpolate and maxr_interpolate>0 and !kspace) return 1.0;
            result = interpolator->Evaluate(r);
        }
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
/*
    if (interpolator_y > 0)
        cerr << "Interpolator was initialized but can't use it at y=" << y << " "
        << LINEINFO << endl;
*/
    /// Initialize new interpolator and use it
    int rind = FindIndex(r, rvals);
    int yind = FindIndex(y, yvals);
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
    double *tmparray = new double[interpo_points];
    double *tmpxarray = new double[interpo_points];
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
    double result=0;
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

bool AmplitudeLib::InterpolatorInitialized(double y)
{
    if (std::abs(y - interpolator_y) < 0.01)
        return true;
    else
        return false;
}

double AmplitudeLib::S(double r, double y, int der)
{
    double s = 1.0 - N(r,y,der);
    if (s<=0) return 0;
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
    double y; double kt;
    bool adjoint;
    AmplitudeLib* N;
};
double N_k_helperf(double r, void* p);
double AmplitudeLib::N_k(double kt, double y)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    bool tmp_range = SetOutOfRangeErrors(false);
    
    N_k_helper par;
    par.y=y; par.N=this; par.kt=kt;
    double result = fourier_j0(kt,N_k_helperf,&par);

    SetOutOfRangeErrors(tmp_range);
    return result;
}

double N_k_helperf(double r, void* p)
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
    double y; double x;
    AmplitudeLib* N;
};
double N_k_to_x_helperf(double k, void* p)
{
    N_k_to_x_helper* par = (N_k_to_x_helper*)p;
    return par->N->N(k, par->y)*k;
}
double AmplitudeLib::N_k_to_x(double x, double y)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    bool tmp_range = SetOutOfRangeErrors(false);
    
    N_k_to_x_helper par;
    par.y=y; par.N=this; par.x=x;
    double result = fourier_j0(x,N_k_to_x_helperf,&par);
    return x*x*result;

    SetOutOfRangeErrors(tmp_range);
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
double S_k_helperf(double r, void* p);
double AmplitudeLib::S_k(double kt, double y, bool adjoint)
{
	SetOutOfRangeErrors(false);
	
	if (!InterpolatorInitialized(y))
		cerr << "Interpolator is not initialized and we are calculating S_k, are you sure? "
			<< "(interpolator y: " << interpolator_y << ", asked y=" << y <<") "
			<< LINEINFO << endl;
	
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    
    N_k_helper par;
    par.y=y; par.N=this; par.kt=kt; par.adjoint=adjoint;

    double result=0;

    if (kt < 1e-3 or true)  // k_T \approx 0 -> integrate just \int d^2 r S(r)
    {
        gsl_function fun; fun.function=S_k_helperf;
        fun.params=&par;
    

        double abserr; size_t eval;
        gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1000);
		int status = gsl_integration_qag(&fun, MinR(), MaxR(), 0, 0.001,
			1000, GSL_INTEG_GAUSS51, ws, &result, &abserr);
		gsl_integration_workspace_free(ws);
        //int status = gsl_integration_qng(&fun, MinR(), MaxR(),
        //    0, 0.01,  &result, &abserr, &eval);
        if (status)
        {
            cerr << "S_k integral failed at " << LINEINFO <<", y=" << y
            << " k_T=" << kt << ", adjoint: " << adjoint << " relerr " << std::abs(abserr/result) << endl;
        }
    }
    else
        result = fourier_j0(kt,S_k_helperf,&par);
    return result*2.0*M_PI; 
}

double S_k_helperf(double r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    double result;
    if (!par->adjoint)
    {
        if (r > par->N->MaxR()) return 0;
        result = r*(1.0-par->N->N(r, par->y));
    }
    else
    {
        result = r*(1.0-par->N->N_A(r, par->y));
    }
    result *= gsl_sf_bessel_J0(par->kt*r);

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
    double Qsqr,y;
    WaveFunction* wf;
};

double Inthelperf_totxs(double r, void* p)
{
    Inthelper_totxs* par = (Inthelper_totxs*)p;

    double result = r*par->N->N(r,par->y);
    if (par->pol==0)    // Longitudinal
        result *= par->wf->PsiSqr_L_intz(par->Qsqr, r);
    else if (par->pol==1)   // Transverse
        result *= par->wf->PsiSqr_T_intz(par->Qsqr, r);
    else
        cerr << "Invalid polarization " << par->pol << " at " << LINEINFO << endl;

    return result;
}

double AmplitudeLib::ProtonPhotonCrossSection(double Qsqr, double y, int pol)
{
    Inthelper_totxs par; par.N=this;
    par.pol=pol; par.Qsqr=Qsqr; par.y=y;
    VirtualPhoton wavef;
    par.wf=&wavef;

    gsl_function fun; fun.function=Inthelperf_totxs;
    fun.params=&par;

    SetOutOfRangeErrors(false);
    

    double result,abserr; size_t eval;
    const int MAXITER_RINT=1000;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_RINT);
    int status = gsl_integration_qag(&fun, 0.1*MinR(), 100.0*MaxR(), 0, 0.001,
        MAXITER_RINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    //int status = gsl_integration_qng(&fun, MinR(), MaxR(),
    //0, 0.001,  &result, &abserr, &eval);
    
    if(status){ std::cerr<< "r integral in ProtonPhotonCrossSection failed with code " 
        << status << " (Qsqr=" << Qsqr << ", y=" << y << " result " << result  
        << " relerr=" << abserr/result << ") at " << LINEINFO << std::endl;
    }

    return 2.0*M_PI*result; //2\pi from \theta integral

    

}


/*
 * d ln N / d ln r^2 = 1/N * d N / d ln r^2 = 1/N d N / dr^2  dr^2 / d ln r^2
 *  = 1/N d N / dr^2 r^2 = 1/N d N / dr  dr / dr^2  r^2 = r/(2N) * dN/dr
 */
double AmplitudeLib::LogLogDerivative(double r, double y)
{
    double dndr = N(r,y,1);
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
    double y;
    double Ns;
    double gammac;
    AmplitudeLib* N;
};
double SatscaleSolverHelperf(double r, void* p)
{
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    return par->N->N(r, par->y) - par->Ns;
}

double SatscaleSolverHelperf_k(double kt, void* p)
{
   
    SatscaleSolverHelper* par = (SatscaleSolverHelper*)p;
    double val = -std::pow(kt, 2.0*par->gammac)
        *par->N->N(kt, par->y);
    return val;
    
}


double AmplitudeLib::SaturationScale(double y, double Ns)
{
    
    SatscaleSolverHelper helper;
    helper.y=y; helper.Ns=Ns; helper.N=this; helper.gammac=0.5;
    const int MAX_ITER = 1000;
    const double ROOTFINDACCURACY = 0.00001;
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


    double sat = gsl_root_fsolver_root(s);

    gsl_root_fsolver_free(s);

    if (!kspace)
        return sat;
    
    f.function=&SatscaleSolverHelperf_k;
  
    const gsl_min_fminimizer_type *Tm = gsl_min_fminimizer_brent;
    gsl_min_fminimizer *sm = gsl_min_fminimizer_alloc(Tm);

    // Interpolation doesn't work well at the edge of the ktsqr range, so limit
    // the studied interval
    double interval_min = MinR()*1.0001;
    double interval_max = MaxR()*0.999;

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
double AmplitudeLib::N_A(double r, double y, int der)
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
    double result=0;

    if (der==0)
    {
        double n = N(r,y);
        result = 2.0*n - n*n;
    }
    else if (der==1)
    {
        double n = N(r,y);
        double dn = N(r,y,1);
        result= 2.0*dn - 2.0*n*dn;
    }
    else if (der==2)
    {
        double n = N(r,y);
        double dn = N(r,y,1);
        double d2n = N(r,y,2);
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

    tmprarray = new double[data.RPoints()];
    tmpnarray = new double[data.RPoints()];
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = std::log(minr*std::pow(rmultiplier, i));
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

    interpolator_y=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)

    //double satscale = 1.0/SaturationScale(0, 0.22);

    /*
    cout << "# Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << " x0 " << X0()
        << " Q_{s,0}^2 = " << SQR(satscale) << " GeV^2 [ N(r=1/Q_s) = 0.22]" << endl;
    */
}

/*
 * Initializes interpolation method with all read data points at given y
 * If bspline is true, then use bspline interpolation, default is false
 */
void AmplitudeLib::InitializeInterpolation(double y, bool bspline)
{
	if (y<0)
	{
		cerr << "Asked to initialize interpolator with negative rapidity! "
			<< "Dont know what to do, panicking... " << LINEINFO << endl;
		exit(1);
	}
	
    if (std::abs(interpolator_y - y) < 0.01) return;    // Already done it
    if (interpolator_y>=0)
    {
        interpolator_y=-1;
        delete interpolator;
    }
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = tmprarray[i];
        if (i==0) tmpr*=1.0001; if (i==rpoints-1) tmpr*=0.9999;
        tmpnarray[i] = N(tmpr, y);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);
    if (bspline)
        interpolator->SetMethod(INTERPOLATE_BSPLINE);
    interpolator->Initialize();
    interpolator_y = y;
    
    maxr_interpolate=-1;
    int iter=0;
    // Find r s.t. N(r)==1
    double step=2; double prevr=0.1;
    const int MAXITER=40;
    for (double r=0.01; r<=MaxR(); r+=step)
    {
        if (N(r,y)>=0.999999)
        {
            if (step<1e-2)
            {
                maxr_interpolate = r;
                break;
            }
            step /= 1.5;
            r=prevr;
        }
        prevr=r;
        iter++;
        if (iter > MAXITER)
        {
            cerr << "Didn't find maxr_interpolate at y=" << y <<", line " << LINEINFO 
                << ", best guess " << r << endl;
            break;  // Didn't find, dont force any upper limit
        }
    }
}



/*
 * Initializes interpolator and returns it
 */
Interpolator* AmplitudeLib::MakeInterpolator(double y)
{
    double* interp_narray = new double[rpoints];
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
double AmplitudeLib::MinR()
{
    return minr;
}

double AmplitudeLib::MaxR()
{
    return minr*std::pow(rmultiplier, rpoints-1);
}

int AmplitudeLib::YVals()
{
    return yvals.size();
}

double AmplitudeLib::MaxY()
{
    return yvals[yvals.size()-1];
}

bool AmplitudeLib::SetOutOfRangeErrors(bool er)
{
    bool outofrange = out_of_range_errors;
    out_of_range_errors=er;
    return outofrange;    
}

double AmplitudeLib::X0()
{
    return x0;
}

