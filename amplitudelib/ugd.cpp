#include "amplitudelib.hpp"
#include "single_inclusive.hpp"
#include "virtual_photon.hpp"
#include "../tools/tools.hpp"
#include "../tools/config.hpp"
#include <string>
#include <cmath>
#include <algorithm>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_sf_bessel.h>
#include "../pdf/pdf.hpp"
#include "../pdf/cteq.hpp"

using namespace Amplitude;

extern "C"
{
    #include "../fourier/fourier.h"
}

SingleInclusive::SingleInclusive(AmplitudeLib* N_)
{
    N=N_;
}


/*
 * k_T factorization result for the dN/(d^2kt dy)
 * \alpha_s/S_T 2/C_f 1/k_T^2 \int d^2 q UGD(q)/q^2 UQD(k-q)/|k-q|^2
 * UGD is AmplitudeLib::UGD()
 * 
 * If N2!=NULL, it is used to compute \phi_2
 */
struct Inthelper_ktfact{ double x1, x2, pt, qt, scale; AmplitudeLib *N1, *N2; };
double Inthelperf_ktfact_q(double q, void* p);
double Inthelperf_ktfact_phi(double phi, void* p);
const int INTPOINTS_KTFACT = 3;	
double SingleInclusive::dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2, double scale )
{
	double x1 = pt*std::exp(-y)/sqrts;
	double x2 = pt*std::exp(y)/sqrts;
	double y1 = std::log(N->X0()/x1);
	double y2;
	if (N2==NULL)
		y2 = std::log(N->X0()/x2);
	else
		y2 = std::log(N2->X0()/x2);
	
	if (y1<0 or y2<0)
	{
		cerr << "Evolution variables y1=" << y1 <<", y2=" << y2 <<" is not possible to compute, too large x. pt=" << pt << ", y=" << y << ", sqrts=" << sqrts <<" " << LINEINFO << endl;
		return 0;
		if (y1<0) y1=0; if (y2<0)y2=0;
		
	}
	
	Inthelper_ktfact par; par.pt=pt; par.N1=N;
    par.x1=x1; par.x2=x2;
	par.N2=N2;
    if (scale<0)
        par.scale=pt*pt;
    else
        par.scale=scale;
	gsl_function fun; fun.params=&par;
	fun.function=Inthelperf_ktfact_q;
	
	if (N2!=NULL) // Initialize interpolators only once
	{
		N->InitializeInterpolation(x1);
		N2->InitializeInterpolation(x2);
	}
	
	double maxq = std::max(4*pt, 40.0);
    if (maxq>80)
        maxq=80;

    
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.001, maxq, 0, 0.05,		
		INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		//#pragma omp critical
		//cerr << "kt-factorization q integral failed at " << LINEINFO <<", pt=" << pt <<", x1=" << x1 <<", x2=" << x2 
		//	<< ", result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    result *= 2.0/(qcd.Cf()*SQR(pt));
    double alphas=qcd.Alphas(scale);
    result *= alphas;
    
    return result;
}

double Inthelperf_ktfact_q(double q, void *p)
{
	Inthelper_ktfact* par = (Inthelper_ktfact*)p;
	par->qt=q;
	
	gsl_function fun; fun.function=Inthelperf_ktfact_phi;
	fun.params=par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(1);
	int status = gsl_integration_qag(&fun, 0, M_PI, 0, 0.05,
		1, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
    result *= 2.0;	// as we integrate [0,pi], not [0, 2pi]   
   /* if (status)
    {
		cerr << "kt-factorization phi integral failed at " << LINEINFO <<", q=" << q
			<< ", k_T=" << par->pt << " relerr " << std::abs(abserr/result) << endl;
    }*/
	return result;
}

double Inthelperf_ktfact_phi(double phi, void* p)
{
	Inthelper_ktfact* par = (Inthelper_ktfact*)p;

	
	double kt_m_qt = std::sqrt( SQR(par->pt) + SQR(par->qt) - 2.0*par->pt*par->qt*std::cos(phi));
	
	// In this limit we would have ugd(0)/0 -> 0
	if (kt_m_qt < 1e-5)
		return 0;	
	
	
	
	double ugd1,ugd2;
	if (par->N2==NULL)
	{
		par->N1->InitializeInterpolation(par->x1);
		ugd1 = par->N1->Dipole_UGD(par->qt, par->x1, par->scale);
		par->N1->InitializeInterpolation(par->x2);
		ugd2 = par->N1->Dipole_UGD(kt_m_qt, par->x2, par->scale);
	} 
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				ugd1 = par->N1->Dipole_UGD(par->qt, par->x1, par->scale);
			}
			#pragma omp section
			{
				ugd2 = par->N2->Dipole_UGD(kt_m_qt, par->x2, par->scale);
			}
		}
		
	}
	return par->qt * ugd1/SQR(par->qt) * ugd2/SQR(kt_m_qt);
	
	
}


/*
 * Single inclusive hadron production in k_T factorization, convoluted
 * with fragmentation function
 * Integrates \int dz/z^2 dN(k_T/z)
 * 
 * Lower cutoff is basically arbitrary now, but in principle it should have
 * little effect in the final result as p_T specturm is steeply falling
 */
 
const int INTPOITNS_KTFACT_Z = 1;
struct Inthelper_ktfact_fragfun{ SingleInclusive* sinc; AmplitudeLib* N2; FragmentationFunction* fragfun; double scale, pt, y, sqrts; Hadron final; PDF *pdf; };
double Inthelperf_ktfact_fragfun(double z, void* p);

double SingleInclusive::dHadronMultiplicity_dyd2pt_ktfact(double y, double pt, double sqrts, FragmentationFunction* fragfun, Hadron final, AmplitudeLib* N2 )
{
	Inthelper_ktfact_fragfun par;
	par.sinc=this; par.y=y; par.y=y; par.pt=pt; par.fragfun=fragfun; par.final=final;
	par.N2=N2;
	par.sqrts=sqrts;
    par.scale=pt*pt;
	
	gsl_function fun; fun.function=Inthelperf_ktfact_fragfun;
	fun.params=&par;
	
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOITNS_KTFACT_Z);
    double minz = std::max(0.05, pt*std::exp(y)/sqrts);
	int status = gsl_integration_qag(&fun, minz, 1.0, 0, 0.05,
		INTPOITNS_KTFACT_Z, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		#pragma omp critical
		cerr << "kt-factorization z integral failed at " << LINEINFO <<", pt=" << pt
			<< " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    return result;
}


double Inthelperf_ktfact_fragfun(double z, void* p)
{
	Inthelper_ktfact_fragfun* par = (Inthelper_ktfact_fragfun*)p;
	double kt = par->pt/z;
	
	double dn = par->sinc->dHadronMultiplicity_dyd2pt_ktfact_parton(par->y, kt, par->sqrts, par->N2, par->scale);

	return 1.0/SQR(z) * dn * par->fragfun->Evaluate(G, par->final, z, std::sqrt(par->scale) );
}
