#include "dis.hpp"
#include "wave_function.hpp"
#include "virtual_photon.hpp"
#include <gsl/gsl_integration.h>

using namespace Amplitude;

DIS::DIS(AmplitudeLib* amp)
{
    N=amp;
}

/*
 * Virtual photon-proton cross sections
 * Separately for transversially and longitudinally polarized photons
 */
struct Inthelper_totxs
{
    AmplitudeLib* N;
    Polarization pol;    // L (longitudinal) or T (transverse)
    double Qsqr,xbj;
    WaveFunction* wf;
};

double Inthelperf_totxs(double r, void* p)
{
    Inthelper_totxs* par = (Inthelper_totxs*)p;

    double result = r*par->N->N(r,par->xbj);
    if (par->pol==L)    // Longitudinal
        result *= par->wf->PsiSqr_L_intz(par->Qsqr, r);
    else if (par->pol==T)   // Transverse
        result *= par->wf->PsiSqr_T_intz(par->Qsqr, r);
    else
        cerr << "Invalid polarization " << par->pol << " at " << LINEINFO << endl;

    return result;
}

/*
 * Compute total gamma-p cross section
 * by default p=LIGHT which means that we sum over light quarks (u,d,s)
 * If p is something else, use the given quark
 *
 * Default value for the mass is -1 which means that the default values
 * for the quark masses (from VirtualPhoton class) are used
 */

double DIS::ProtonPhotonCrossSection(double Qsqr, double xbj, Polarization pol,Parton p, double mass)
{
    Inthelper_totxs par; par.N=N;
    par.pol=pol; par.Qsqr=Qsqr; par.xbj=xbj;

    
    VirtualPhoton wavef;

    wavef.SetQuark(p, mass);
    
    par.wf=&wavef;

    gsl_function fun; fun.function=Inthelperf_totxs;
    fun.params=&par;

    N->SetOutOfRangeErrors(false);
    

    double result,abserr; 
    const int MAXITER_RINT=100;
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(MAXITER_RINT);
    int status = gsl_integration_qag(&fun, 0.1*N->MinR(), 10.0*N->MaxR(), 0, 0.01,
        MAXITER_RINT, GSL_INTEG_GAUSS51, ws, &result, &abserr);
    gsl_integration_workspace_free(ws);
    //int status = gsl_integration_qng(&fun, MinR(), MaxR(),
    //0, 0.001,  &result, &abserr, &eval);
    
    if(status){ std::cerr<< "r integral in ProtonPhotonCrossSection failed with code " 
        << status << " (Qsqr=" << Qsqr << ", xbj=" << xbj << " result " << result  
        << " relerr=" << abserr/result << ") at " << LINEINFO << std::endl;
    }

    return 2.0*M_PI*result; //2\pi from \theta integral

}

double DIS::F2(double qsqr, double xbj, Parton p, double mass)
{
	double xs_l, xs_t;
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			xs_l = ProtonPhotonCrossSection(qsqr, xbj, L, p, mass);
		}
		#pragma omp section
		{
			xs_t = ProtonPhotonCrossSection(qsqr, xbj, T, p, mass);
		}
	}
	
    return qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
}

double DIS::FL(double qsqr, double xbj, Parton p, double mass)
{
	return qsqr/(4.0*SQR(M_PI)*ALPHA_e) * ProtonPhotonCrossSection(qsqr, xbj, L, p, mass);	
}

double DIS::ReducedCrossSection(double qsqr, double xbj, double sqrts, Parton p, double mass)
{
	double kin_y = qsqr/(sqrts*sqrts*xbj);   // inelasticity, not rapidity
	
	if (xbj > N->X0())
	{
		cerr << "Asked to calculate DIS at x=" << xbj << ", Q^2=" << qsqr << ", mass=" << mass << ", evaluating at x=x0=" << N->X0() << endl;
		xbj = N->X0();
	}
	N->InitializeInterpolation(xbj);
	
	
	double xs_l, xs_t;
	#pragma omp parallel sections
	{
		#pragma omp section
		{
			xs_l = ProtonPhotonCrossSection(qsqr, xbj, L, p, mass);
		}
		#pragma omp section
		{
			xs_t = ProtonPhotonCrossSection(qsqr, xbj, T, p, mass);
		}
	}
	double f2 = qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t);
	double fl = qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l;
		
	return f2 - SQR(kin_y) / ( 1.0 + SQR(1.0-kin_y) ) * fl;
}
