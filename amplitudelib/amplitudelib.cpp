/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
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
#include <string>
#include <sstream>


#include <algorithm>

using namespace Amplitude;


extern "C"
{
    #include "../fourier/fourier.h"
}

/*
 * Load data from a given file
 * Format is specified in file bk/README
 * If kspace is true, the datafile is in kspace solved using BK code
 * ("x axis" in kt^2 in datafile, don't limit N<=1)
 */
AmplitudeLib::AmplitudeLib(std::string datafile, bool kspace_)
{
    // Read BK solution
    DataFile data(datafile);
    data.GetData(n, yvals);

    // Initialize
    out_of_range_errors=true;
    kspace=kspace_;
    
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

    interpolator_xbj=-1;  // if >=0, interpolator is initialized, must free
    // memory (delete tmprarray and tmpnarray at the end)


	std::stringstream ss;
    ss << "Data read from file " << datafile << ", minr: " << minr
        << " maxr: " << MaxR() << " rpoints: " << rpoints << " maxy "
        << yvals[yvals.size()-1] << " x0 " << X0()
        << " Q_{s,0}^2 = " << 2.0/SQR(SaturationScale(x0, 0.393469)) << " GeV^2 [ N(r^2=2/Q_s^2, x=x0) = 0.3934]"
        << " (AmplitudeLib v. " << AMPLITUDELIB_VERSION << ")" ;
    info_string = ss.str();
}

/*
 * Release reserved memory
 */

AmplitudeLib::~AmplitudeLib()
{
    if (interpolator_xbj>=0)
    {
        delete interpolator;
    }
    delete[] tmpnarray;
    delete[] tmprarray;
}

/*
 * Calculate amplitude interpolating data
 * Interpolate rapidity linearly and r using spline
 */
double AmplitudeLib::N(double r, double xbj)
{    
    if (r < MinR() or r > MaxR() )
    {
        if (out_of_range_errors)
            cerr << "r must be between limits [" << MinR() << ", " << MaxR() << "]"
                << " asked r=" << r << " " << LINEINFO
                << endl;
        if (r<MinR()) r=MinR()*1.000001; else if (r>MaxR()) r=MaxR()*0.999999;
    }

    double y = std::log(X0()/xbj);
    
    if (y<0 or y>yvals[yvals.size()-1] )
    {
        if (out_of_range_errors)
            cerr << "y must be between limits [" << 0 << ", "
                << yvals[yvals.size()-1] << "], asked y=" << y << " "
                << LINEINFO << endl;
        if (y < 0) y=0; else if (y>yvals[yvals.size()-1]) y=yvals[yvals.size()-1];
    }

    
    // Use already initialized interpolator
    if (std::abs(xbj - interpolator_xbj)/xbj < 0.001 )
    {
        double result=0;

        // Can't interpolate (too large dipole), return 1.0
        if (r >= maxr_interpolate and maxr_interpolate>0 and !kspace) { return 1.0; }

        result = interpolator->Evaluate(r);
        
        if (!kspace and result>1) return 1;  // in x space limit N(r)<=1
        if (result<0) return 0;              // limit N >= 0
        return result;

    }

    /// Initialize new interpolator and use it
    int rind = FindIndex(r, rvals);
    int yind = FindIndex(y, yvals);
    int interpolation_points = INTERPOLATION_POINTS;

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
    interp.Initialize();
    double result=0;

    result = interp.Evaluate(r);

    
    delete[] tmparray;
    delete[] tmpxarray;

    return result;

}

bool AmplitudeLib::InterpolatorInitialized(double xbj)
{
    if (std::abs(xbj - interpolator_xbj)/xbj < 0.001)
        return true;
    else
        return false;
}

double AmplitudeLib::S(double r, double xbj)
{
    double s = 1.0 - N(r,xbj);
    if (s<=0) return 0;
    if (s>=1.0) return 1.0;
    return s;
}

/*
 * FT the amplitude to the k-space (WW UGD)
 * N(k) = \int d^2 r/(2\pi r^2) exp(ik.r)N(r)
 *  = \int dr/r BesselJ[0,k*r] * N(r)
 * 
 * Note: for performance reasons it is probably a good idea to
 * call AmplitudeLib::InitializeInterpolation(y) before this
 */
struct N_k_helper
{
    double xbj; double kt; double power;
    bool adjoint;
    AmplitudeLib* N;
};
double N_k_helperf(double r, void* p);
double AmplitudeLib::WW_UGD(double kt, double xbj)
{
    // Some initialisation stuff
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    bool tmp_range = SetOutOfRangeErrors(false);
    
    N_k_helper par;
    par.xbj=xbj; par.N=this; par.kt=kt;
    double result = fourier_j0(kt,N_k_helperf,&par);
    /// \todo Respect FOURIER_TRANSFER parameter
    SetOutOfRangeErrors(tmp_range);
    return result;
}

double N_k_helperf(double r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    if (r < par->N->MinR()) return 0;
    else if (r > par->N->MaxR()) return 1.0/r;
    return 1.0/r*par->N->N(r, par->xbj);
}

/*
 * FT amplitude from k to x space
 * N(r) = r^2 \int d^2 k/(2\pi) exp(-ik.r) N(k)
 *      = r^2 \int dk k BesselJ[0,k*r] * N(k)
 */
struct N_k_to_x_helper
{
    double xbj; double x;
    AmplitudeLib* N;
};
double N_k_to_x_helperf(double k, void* p)
{
    N_k_to_x_helper* par = (N_k_to_x_helper*)p;
    return par->N->N(k, par->xbj)*k;
}
double AmplitudeLib::N_k_to_x(double r, double xbj)
{
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    bool tmp_range = SetOutOfRangeErrors(false);
    
    N_k_to_x_helper par;
    par.xbj=xbj; par.N=this; par.x=r;
    double result = fourier_j0(r,N_k_to_x_helperf,&par);
    return r*r*result;

    SetOutOfRangeErrors(tmp_range);
}


/*
 * FT S=1-N to the k-space
 * S(k) = \int d^2 r exp(ik.r) (1-N(r))^power
 *  = (2\pi) \int dr r BesselJ[0,k*r] * (1-N(r))^power
 * 
 * Note: for performance reasons it is probably a good idea to
 * call AmplitudeLib::InitializeInterpolation(y) before this
 *
 * if adjoint==true, then use adjoint representation for N
 * Default value for adjoint is false
 * 
 * Default value for power is 1.0
 */


double S_k_helperf(double r, void* p);
double AmplitudeLib::S_k(double kt, double xbj, Representation rep, double pow)
{
	SetOutOfRangeErrors(false);
	if (!InterpolatorInitialized(xbj))
		cerr << "Interpolator is not initialized and we are calculating S_k, are you sure? "
			<< "(interpolator y: " << interpolator_xbj << ", asked x=" << xbj <<") "
			<< LINEINFO << endl;
	
    // Some initialisation stuff -
    set_fpu_state();
    init_workspace_fourier(FOURIER_ZEROS);   // number of bessel zeroes, max 2000
    
    N_k_helper par;
    par.xbj=xbj; par.N=this; par.kt=kt;  par.power=pow;

    if (rep==ADJOINT)
        par.adjoint=true;
    else
        par.adjoint=false;

    double result=0;

    if (kt < 1e-3  or GetFTMethod()==GSL)  // k_T \approx 0 -> integrate just \int d^2 r S(r)^power
    {
        gsl_function fun; fun.function=S_k_helperf;
        fun.params=&par;
    

        double abserr; 
        gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1000);
		int status = gsl_integration_qag(&fun, MinR(), MaxR(), 0, 0.001,
			1000, GSL_INTEG_GAUSS51, ws, &result, &abserr);
		gsl_integration_workspace_free(ws);
        //int status = gsl_integration_qng(&fun, MinR(), MaxR(),
        //    0, 0.01,  &result, &abserr, &eval);
        if (status)
        {
            cerr << "S_k integral failed at " << LINEINFO <<", x=" << xbj
            << " k_T=" << kt << ", adjoint: " << par.adjoint << " relerr " << std::abs(abserr/result) << endl;
        }
    }
    else
        result = fourier_j0(kt,S_k_helperf,&par);
        
        
    if (result<-0.00001){
		cerr << "S_k transfer is negative=" << result*2.0*M_PI <<", kt=" << kt <<", x=" << xbj <<", adj: " << par.adjoint <<", power=" << pow << " " << LINEINFO << endl;
    }
    return result*2.0*M_PI; 
}

double S_k_helperf(double r, void* p)
{
    N_k_helper* par = (N_k_helper*) p;
    double result;
    if (!par->adjoint)
    {
        if (r > par->N->MaxR()) return 0;
        result = r*std::pow(par->N->S(r, par->xbj), par->power);
    }
    else
    {
        result = r*std::pow(1.0-par->N->N_A(r, par->xbj), par->power);
    }
    
    // J0 is in fourier_j0() fun
    if (par->N->GetFTMethod()==GSL)		// not using fourier_j0
		result *= gsl_sf_bessel_J0(par->kt*r);

    return result;
}




/*
 * d ln N / d ln r^2 = 1/N * d N / d ln r^2 = 1/N d N / dr^2  dr^2 / d ln r^2
 *  = 1/N d N / dr^2 r^2 = 1/N d N / dr  dr / dr^2  r^2 = r/(2N) * dN/dr
 */
double AmplitudeLib::LogLogDerivative(double r, double y)
{
    ///TODO: Requires derivative
    /*double dndr = N(r,y,1);
    return r/(2.0*N(r,y))*dndr;*/
    return 0;
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

    if (gsl_min_fminimizer_set(sm, &f, sat, interval_min, interval_max)
        == GSL_EINVAL)
    {
        // Endpoints do not enclose a minimum (on the other hand
        // helper(pos) > helper(min),helper(max), but we can anyway continue
        cerr << "Endpoints do not enclose a minumum! " << LINEINFO << endl;
    }


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
double AmplitudeLib::N_A(double r, double y)
{
    if (kspace)
    {
        cerr <<"N_A is not implemented for mome. space!" << endl;
        return 0;
    }

    double result=0;

    double n = N(r,y);
    result = 2.0*n - n*n;
    if (result>1) result=1;
    if (result<0) result=0;
    return result;
    
    

    return result;
}




/*
 * Initializes interpolation method with all read data points at given y
 * If bspline is true, then use bspline interpolation, default is false
 */
void AmplitudeLib::InitializeInterpolation(double xbj)
{
	if (xbj > X0())
	{
		cerr << "Asked to initialize interpolator with too large x=" << xbj
        << " (x0=" << X0() <<")! "
        << "Dont know what to do, panicking... " << LINEINFO << endl;
		exit(1);
	}

    double y = std::log(xbj/X0());
	if (y>MaxY())
	{
		cerr << "Asked to initialize interpolator with too large y=" << y <<", maxy=" << MaxY() << endl;
		exit(1);
	}
	
    if (std::abs((interpolator_xbj - xbj)/xbj) < 0.001) return;    // Already done it
    if (interpolator_xbj>=0)
    {
        interpolator_xbj=-1;
        delete interpolator;
    }
    for (int i=0; i<rpoints; i++)
    {
        double tmpr = tmprarray[i];
        if (i==0) tmpr*=1.0001; if (i==rpoints-1) tmpr*=0.9999;
        tmpnarray[i] = N(tmpr, xbj);
    }
    interpolator = new Interpolator(tmprarray, tmpnarray, rpoints);

    interpolator->Initialize();
    interpolator_xbj = xbj;
    
    maxr_interpolate=-1;
    int iter=0;
    // Find r s.t. N(r)==1
    double step=2; double prevr=0.1;
    const int MAXITER=40;
    for (double r=0.01; r<=MaxR(); r+=step)
    {
        if (N(r,y)>=0.99999)		// check that this accuracy is the same as in rbk/src/solver.cpp
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
            //cerr << "Didn't find maxr_interpolate at y=" << y <<", ignoring it " << LINEINFO << endl;
            maxr_interpolate=-1;
            break;  // Didn't find, dont force any upper limit
        }
    }
}



/*
 * Initializes interpolator and returns it
 */
Interpolator* AmplitudeLib::MakeInterpolator(double xbj)
{
    std::vector<double> tmpnvals;
    for (int i=0; i<rpoints; i++)
    {
        tmpnvals.push_back(N(rvals[i], xbj));
    }
    Interpolator* inter = new Interpolator(rvals, tmpnvals);
    inter->Initialize();

    return inter;
        
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

int AmplitudeLib::YPoints()
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

void AmplitudeLib::SetX0(double x0_)
{
	x0=x0_;
}

std::string AmplitudeLib::GetString()
{
	return info_string;
}

void AmplitudeLib::SetFTMethod(FT_Method f)
{
	ft=f;	
}


FT_Method AmplitudeLib::GetFTMethod()
{
	return ft;
}

std::string AmplitudeLib::Version()
{
	std::stringstream s;
	s << AMPLITUDELIB_VERSION << " (build " <<  __DATE__ << " " << __TIME__ << ")";
	return s.str();
}

/*
 * UGD normalization here is the "KMR" normalization, see e.g. 
 * Ref. hep-ph/0101348] (eq (26) and hep-ph/0111362 (eq (41)
 */


/* Dipole (KMR) UGD C_F/(8\pi^3) S_T/\alpha_s(q) q^4 S_k(q)
 * default value of S_T is 1.0, so it is left for the user to 
 * specify
 * Scale is the scale at which alpha_s is evaluated, if negative (default)
 * use q^2
 * S_k is 2d FT of S(r) without extra 2pi factors,
 * S(k) = \int d^2 r e^(iqr) S(r), AmplitudeLib::S_
 */
double AmplitudeLib::Dipole_UGD(double q, double xbj, double scale_, double S_T)
{
	double scale;
	if (scale_<0) scale=q*q; else scale=scale_;
	double alphas = qcd.Alphas(scale);
	return qcd.Cf() / (8.0 * M_PI*M_PI*M_PI) * S_T/alphas * std::pow(q,4) * S_k(q, xbj, ADJOINT);

}

/*
 * Integrated GD, x*g(x), from UGD
 * The Q^2 dependence comes as a upper limit of an integral
 * xg(x,Q^2) = \int_0^Q dq^2/q^2 UGD(q)
 */
double Inthelperf_xg(double qsqr, void* p);
struct Inthelper_xg
{
	double xbj,q; AmplitudeLib* N;
};
double AmplitudeLib::xg(double x, double q)
{
	Inthelper_xg par; 
	par.N=this;
    InitializeInterpolation(x);
	par.xbj=x; par.q=q;
	
	gsl_function fun; fun.function=Inthelperf_xg;
	fun.params=&par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(100);
	int status = gsl_integration_qag(&fun, 0, SQR(q), 0, 0.01,
		100, GSL_INTEG_GAUSS51, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		cerr << "UGD integral failed at " << LINEINFO <<", x=" << x
			<< ", k_T=" << q << " res " << result << " relerr " << std::abs(abserr/result)  << endl;
    }
	return result;
}

double Inthelperf_xg(double qsqr, void* p)
{
	Inthelper_xg* par = (Inthelper_xg*) p;
	return 1.0/qsqr * par->N->Dipole_UGD(std::sqrt(qsqr), par->xbj, SQR(par->q));
}

QCD& AmplitudeLib::GetQCD()
{
    return qcd;
}

