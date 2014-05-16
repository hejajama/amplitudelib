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
#include "../pdf/pdf.hpp"
#include "../pdf/cteq.hpp"

using namespace Amplitude;

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
	double alphas = Alphas(scale);
	return Cf / (8.0 * M_PI*M_PI*M_PI) * S_T/alphas * std::pow(q,4) * S_k(q, y, true);

}

/*
 * Integrated GD, x*g(x), from UGD
 * The Q^2 dependence comes as a upper limit of an integral
 * xg(x,Q^2) = \int_0^Q dq^2/q^2 UGD(q)
 */
double Inthelperf_xg(double qsqr, void* p);
struct Inthelper_xg
{
	double y,q; AmplitudeLib* N;
};
double AmplitudeLib::xg(double x, double q)
{
	Inthelper_xg par; 
	par.N=this;
	double y = std::log(X0()/x); InitializeInterpolation(y);
	par.y=y; par.q=q;
	
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
	return 1.0/qsqr * par->N->UGD(std::sqrt(qsqr), par->y, SQR(par->q));
}


/*
 * k_T factorization result for the dN/(d^2kt dy)
 * \alpha_s/S_T 2/C_f 1/k_T^2 \int d^2 q UGD(q)/q^2 UQD(k-q)/|k-q|^2
 * UGD is AmplitudeLib::UGD()
 * 
 * If N2!=NULL, it is used to compute \phi_2
 *
 * If scale<0 (default), use parton pt as a alphas scale, otherwise use given scale
 */
struct Inthelper_ktfact{ double y1, y2, pt, qt, scale; AmplitudeLib *N1, *N2; };
double Inthelperf_ktfact_q(double q, void* p);
double Inthelperf_ktfact_phi(double phi, void* p);
const int INTPOINTS_KTFACT = 3;	///DEBUG orig 2
double AmplitudeLib::dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2, double scale )
{
	double x1 = pt*std::exp(-y)/sqrts;
	double x2 = pt*std::exp(y)/sqrts;
	double y1 = std::log(X0()/x1);
	double y2;
	if (N2==NULL)
		y2 = std::log(X0()/x2);
	else
		y2 = std::log(N2->X0()/x2);
	
	if (y1<0 or y2<0)
	{
		#pragma omp critical
		cerr << "y1=" << y1 <<", y2=" << y2 <<" is not possible to compute, too large x. pt=" << pt << ", y=" << y << ", sqrts=" << sqrts <<" " << LINEINFO << endl;
		return 0;
		if (y1<0) y1=0; if (y2<0)y2=0;
		
	}
	
	Inthelper_ktfact par; par.y1=y1; par.y2=y2; par.pt=pt; par.N1=this;
	par.N2=N2;
    if (scale<0)
        par.scale = pt*pt;
    else
        par.scale=scale;
	gsl_function fun; fun.params=&par;
	fun.function=Inthelperf_ktfact_q;
	
	if (N2!=NULL) // Initialize interpolators only once
	{
		InitializeInterpolation(y1);
		N2->InitializeInterpolation(y2);
	}
	
	double maxq = std::max(4*pt, 40.0);
    if (maxq>60)
        maxq=60;
	
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
    
    result *= 2.0/(Cf*SQR(pt));
    double alphas=Alphas(scale);
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
		par->N1->InitializeInterpolation(par->y1);
		ugd1 = par->N1->UGD(par->qt, par->y1, par->scale);
		par->N1->InitializeInterpolation(par->y2);
		ugd2 = par->N1->UGD(kt_m_qt, par->y2, par->scale);
	} 
	else
	{
		#pragma omp parallel sections
		{
			#pragma omp section
			{
				ugd1 = par->N1->UGD(par->qt, par->y1, par->scale);
			}
			#pragma omp section
			{
				ugd2 = par->N2->UGD(kt_m_qt, par->y2, par->scale);
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
 
 // if true, change between ktfact and hybrid formalism when x_p>X0()
const bool HADRONPROD_TRANSITION = false;
const int INTPOITNS_KTFACT_Z = 1;
struct Inthelper_ktfact_fragfun{ AmplitudeLib* N1, *N2; FragmentationFunction* fragfun; double scale, pt, y, sqrts; Hadron final; PDF *pdf; };
double Inthelperf_ktfact_fragfun(double z, void* p);

double AmplitudeLib::dHadronMultiplicity_dyd2pt_ktfact(double y, double pt, double sqrts, FragmentationFunction* fragfun, Hadron final, AmplitudeLib* N2)
{
    
	Inthelper_ktfact_fragfun par;
	par.N1=this; par.y=y; par.y=y; par.pt=pt; par.fragfun=fragfun; par.final=final;
	par.N2=N2;
    par.scale=std::max(1.0,pt*pt);
	par.sqrts=sqrts;
	
	PDF* pdf;
	if (HADRONPROD_TRANSITION)
	{
		pdf = new CTEQ(); pdf->SetOrder(LO); pdf->Initialize();
		par.pdf=pdf;
	}
	
	gsl_function fun; fun.function=Inthelperf_ktfact_fragfun;
	fun.params=&par;
	
	double result, abserr; 
    gsl_integration_workspace* ws = gsl_integration_workspace_alloc(INTPOITNS_KTFACT_Z);
	int status = gsl_integration_qag(&fun, std::max(0.1,pt*std::exp(y)/sqrts), 1.0, 0, 0.05,
		INTPOITNS_KTFACT_Z, GSL_INTEG_GAUSS15, ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
       
    if (status)
    {
		#pragma omp critical
		cerr << "kt-factorization z integral failed at " << LINEINFO <<", pt=" << pt
			<< " result " << result << " relerr " << std::abs(abserr/result) << endl;
    }
    
    if (HADRONPROD_TRANSITION) delete pdf;
    
    return result;
}


double Inthelperf_ktfact_fragfun(double z, void* p)
{
	Inthelper_ktfact_fragfun* par = (Inthelper_ktfact_fragfun*)p;
	double kt = par->pt/z;
	double scale = std::sqrt(par->scale);	
	double xp = kt * std::exp(par->y) / par->sqrts;
	if (HADRONPROD_TRANSITION and xp > par->N1->X0())
	{
		// Entering the region where hybrid formalism is ok,
		// quite ugly copypaste from xs.cpp:Inthelperf_hadronprod
		cout << "# switch to hybrid, xp=" << xp <<", xa=" << xp*std::exp(-2.0*par->y) <<" kt=" << kt << endl;
		double hybridres=0;
		
		double y_A = std::log(par->N1->X0() / (kt*std::exp(-par->y)/par->sqrts));
		if (y_A<0) return 0;	// Most likely very large pt -> tiny contribution
		par->N1->InitializeInterpolation(y_A);
		double nf = par->N1->S_k(kt, y_A);
		// PDF and fragmentation
		double xqf = par->pdf->xq(xp, scale, U)*par->fragfun->Evaluate(U, par->final, z, scale)
			+ par->pdf->xq(xp, scale, D)*par->fragfun->Evaluate(D, par->final, z, scale)
			+ par->pdf->xq(xp, scale, S)*par->fragfun->Evaluate(S, par->final, z, scale);

		hybridres = nf*xqf;
		
		// Adjoint representation, gluon scatters
		double na = par->N1->S_k(kt, y_A, true);
		double xgf = par->pdf->xq(xp, scale, G)*par->fragfun->Evaluate(G, par->final, z, scale);
		hybridres += na*xgf;
		
		return hybridres/(SQR(2.0*M_PI)*SQR(z));
	}
	double dn = par->N1->dHadronMultiplicity_dyd2pt_ktfact_parton(par->y, kt, par->sqrts, par->N2, par->scale);

	return 1.0/SQR(z) * dn * par->fragfun->Evaluate(G, par->final, z, scale );
}
