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
 * UGD normalization here is the "KMR" normalization, see e.g. 
 * Ref. hep-ph/0101348] (eq (26) and hep-ph/0111362 (eq (41)
 */


/* KMR UGD C_F/(8\pi^3) S_T/\alpha_s(q) q^4 S_k(q)
 * default value of S_T is 1.0, so it is left for the user to 
 * specify
 * Scale is the scale at which alpha_s is evaluated, if negative (default)
 * use q^2
 * S_k is 2d FT of S(r) without extra 2pi factors,
 * S(k) = \int d^2 r e^(iqr) S(r), AmplitudeLib::S_
 */
double AmplitudeLib::UGD(double q, double y, double scale_, double S_T)
{
	double scale;
	if (scale_<0) scale=q*q; else scale=scale_;
	return Cf / (8.0 * M_PI*M_PI*M_PI) * S_T/Alpha_s(scale) * std::pow(q,4) * S_k(q, y, true);

}

/*
 * Integrated GD, x*g(x), from UGD
 * The Q^2 dependence comes as a upper limit of an integral
 * xg(x,Q^2) = \int_0^Q dq^2/q^2 UGD(q)
 */
double Inthelperf_xg(double qsqr, void* p);
struct Inthelper_xg
{
	double y; AmplitudeLib* N;
};
double AmplitudeLib::xg(double x, double q)
{
	Inthelper_xg par; 
	par.N=this;
	double y = std::log(X0()/x); InitializeInterpolation(y);
	par.y=y;
	
	gsl_function fun; fun.function=Inthelperf_xg;
	fun.params=&par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(100);
	int status = gsl_integration_qag(&fun, 0, SQR(q), 0, 0.001,
		100, GSL_INTEG_GAUSS51, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		cerr << "UGD integral failed at " << LINEINFO <<", x=" << x
			<< ", k_T=" << q << " relerr " << std::abs(abserr/result) << endl;
    }
	return result;
}

double Inthelperf_xg(double qsqr, void* p)
{
	Inthelper_xg* par = (Inthelper_xg*) p;
	return 1.0/qsqr * par->N->UGD(std::sqrt(qsqr), par->y);
}


/*
 * k_T factorization result for the dN/(d^2kt dy)
 * \alpha_s/S_T 2/C_f 1/k_T^2 \int d^2 q UGD(q)/q^2 UQD(k-q)/|k-q|^2
 * UGD is AmplitudeLib::UGD()
 * 
 * If N2!=NULL, it is used to compute \phi_2
 */
struct Inthelper_ktfact{ double y1, y2, pt, qt; AmplitudeLib *N1, *N2; };
double Inthelperf_ktfact_q(double q, void* p);
double Inthelperf_ktfact_phi(double phi, void* p);
const int INTPOINTS_KTFACT = 1;
double AmplitudeLib::dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2 )
{
	double x1 = pt*std::exp(y)/sqrts;
	double x2 = pt*std::exp(-y)/sqrts;
	double y1 = std::log(X0()/x1);
	double y2 = std::log(X0()/x2);
	
	if (y1<0 or y2<0)
	{
		cerr << "y1=" << y1 <<", y2=" << y2 <<" is not possible to compute, too small x. pt=" << pt << ", y=" << y << " " << LINEINFO << endl;
		if (y1<0) y1=0; if (y2<0)y2=0;
		
	}
	
	Inthelper_ktfact par; par.y1=y1; par.y2=y2; par.pt=pt; par.N1=this;
	par.N2=N2;
	gsl_function fun; fun.params=&par;
	fun.function=Inthelperf_ktfact_q;
	
	if (N2!=NULL) // Initialize interpolators only once
	{
		InitializeInterpolation(y1);
		N2->InitializeInterpolation(y2);
	}
	
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(2*INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.01, std::max(25.0,7.0*pt), 0, 0.05,
		2*INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		cerr << "kt-factorization q integral failed at " << LINEINFO <<", pt=" << pt
			<< " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    result *= 2.0/(Cf*SQR(pt));
    result *= Alpha_s(pt*pt);
    
    return result;
}

double Inthelperf_ktfact_q(double q, void *p)
{
	Inthelper_ktfact* par = (Inthelper_ktfact*)p;
	par->qt=q;
	
	gsl_function fun; fun.function=Inthelperf_ktfact_phi;
	fun.params=par;
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0, 2.0*M_PI, 0, 0.05,
		INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
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
	
	// Int this limit we would have ugd(0)/0 -> 0
	if (kt_m_qt < 1e-5)
		return 0;	
	
	double ugd1,ugd2;
	if (par->N2==NULL)
	{
		par->N1->InitializeInterpolation(par->y1);
		ugd1 = par->N1->UGD(par->qt, par->y1, SQR(par->pt));
		par->N1->InitializeInterpolation(par->y2);
		ugd2 = par->N1->UGD(kt_m_qt, par->y2, SQR(par->pt));
	} 
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				ugd1 = par->N1->UGD(par->qt, par->y1, SQR(par->pt));
			}
			#pragma omp section
			{
				ugd2 = par->N2->UGD(kt_m_qt, par->y2, SQR(par->pt));
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
struct Inthelper_ktfact_fragfun{ AmplitudeLib* N1, *N2; FragmentationFunction* fragfun; double scale, pt, y, sqrts; Hadron final; };
double Inthelperf_ktfact_fragfun(double z, void* p);
double AmplitudeLib::dHadronMultiplicity_dyd2pt_ktfact(double y, double pt, double sqrts, FragmentationFunction* fragfun, Hadron final, AmplitudeLib* N2 )
{
	Inthelper_ktfact_fragfun par;
	par.N1=this; par.y=y; par.y=y; par.pt=pt; par.fragfun=fragfun; par.final=final;
	par.N2=N2;
	par.sqrts=sqrts;
	
	gsl_function fun; fun.function=Inthelperf_ktfact_fragfun;
	fun.params=&par;
	
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(2*INTPOINTS_KTFACT);
	int status = gsl_integration_qag(&fun, 0.1, 1.0, 0, 0.05,
		2*INTPOINTS_KTFACT, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		cerr << "kt-factorization z integral failed at " << LINEINFO <<", pt=" << pt
			<< " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    return result;
}

double Inthelperf_ktfact_fragfun(double z, void* p)
{
	Inthelper_ktfact_fragfun* par = (Inthelper_ktfact_fragfun*)p;
	double kt = par->pt/z;
	double dn = par->N1->dHadronMultiplicity_dyd2pt_ktfact_parton(par->y, kt, par->sqrts, par->N2);
	return 1.0/SQR(z) * dn * par->fragfun->Evaluate(G, par->final, z, std::max(1.0, par->pt) );
}
