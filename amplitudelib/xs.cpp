/*
 * AmplitudeLib cross section calculation methods
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

#include "amplitudelib.hpp"
#include "../tools/tools.hpp"
#include <cmath>
#include "../pdf/pdf.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../fragmentation/kkp.hpp"
#include "../fragmentation/pkhff.hpp"
#include "../pdf/cteq.hpp"
#include <gsl/gsl_integration.h>

/* Differential forward hadron production multiplicity
 * dN_h / (dy_h d^2 p_T)
 * Ref 1001.1378 eq (1), normalization is arbitrary
 */
struct Inthelper_hadronprod
{   
    AmplitudeLib *N;
    double y, pt, xf;
    double sqrts;
    double miny,maxy;
    double minpt,maxpt;
    PDF* pdf;
    FragmentationFunction* frag;
    bool deuteron;
    Hadron final;
    bool ptint;     // if false, don't integrate over pt
};

double Inthelperf_hadronprod(double z, void *p)
{
    // x1: fraction of parent hadron momentum carried by quark,
    // x2: Bjorken x for the dense system
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    double x1 = par->xf/z; double x2 = x1*std::exp(-2.0*par->y);
    if (x1>1)
    {
		cerr << "Entering kinematically fobidden region at y=" << par->y <<", pt=" << par->pt << " " << LINEINFO << endl;
		return 0;
	}
    
    double y_A = std::log(par->N->X0()/x2);

    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", z " << z << " xf " <<
            par->xf << " x1 " << x1 << " x2 " << x2 << " y_A " << y_A 
            << " y " << par->y << " sqrts " 
            << par->sqrts << " pt " << par->pt << endl ;
        return 0;
    }

    bool deuteron = par->deuteron;
    par->N->InitializeInterpolation(y_A);

    double result = 0;

    // Quark from proton:
    // UGD in fundamental representation
    double nf = par->N->S_k(par->pt/z, y_A);
    // PDF and fragmentation
    double xqf = par->pdf->xq(x1, par->pt, U)*par->frag->Evaluate(U, par->final, z, par->pt)
        + par->pdf->xq(x1, par->pt, D)*par->frag->Evaluate(D, par->final, z, par->pt);

    if (deuteron)
    {
        // isospin symmetry, u in p -> d in n
        xqf += par->pdf->xq(x1, par->pt, U)*par->frag->Evaluate(D, par->final, z, par->pt)
        + par->pdf->xq(x1, par->pt, D)*par->frag->Evaluate(U, par->final, z, par->pt);
    }
        
    result = nf*xqf;

    // Adjoint representation, gluon scatters
    double na = par->N->S_k(par->pt/z, y_A, true);
    xqf = par->pdf->xq(x1, par->pt, G)*par->frag->Evaluate(G, par->final, z, par->pt);
    if (deuteron) xqf *= 2.0;   // gluon pdf gets multiplied by 2
    result += na*xqf;

    return result/SQR(z);
}

/*
 * Calculate dN / dyd^2 pt
 * if deuteron is true (default: false), the probe is deuteron, not proton
 * => use isospin symmetry
 * Default final state particle is PI0
 */
double AmplitudeLib::dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
    FragmentationFunction* fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    const double K = 1.0; // normalization factor
    // We assume light hadrons
    double xf = pt/sqrts*std::exp(y);
    
    SetOutOfRangeErrors(false);
    
    if (xf > 1 or sqrts < 10)
    {
        cerr << "Parameters don't make sense, xf=" << xf << ", sqrts="
            << sqrts << ", y=" << y << " pt=" << pt << " " << LINEINFO << endl;
        return 0;
    }

    Inthelper_hadronprod helper;
    helper.N=this; helper.y=y; helper.pt=pt; helper.xf=xf;
    helper.deuteron=deuteron;
    helper.final=final;
    helper.sqrts=sqrts;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.y=y;

    double result=0; double abserr=0;
    const int MULTIPLICITYXINTPOINTS=5;

    
    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_hadronprod;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MULTIPLICITYXINTPOINTS);

    int status;
    status=gsl_integration_qag(&fun, xf, 1.0,
            0, 0.05, MULTIPLICITYXINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "z integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << endl;


    result *= K / SQR(2.0*M_PI);

    return result;
}

/*
 * Multiplicity integrated over rapidity and p_T range
 */
double Inthelperf_hadronprod_yint(double y, void *p);
double Inthelperf_hadronprod_ptint(double pt, void *p);
const int HADRONPROD_YINTPOINTS=3;
const int HADRONPROD_PTINTPOINTS=3;
const double HADRONPROD_INTACCURACY=0.01;

double AmplitudeLib::HadronMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    Inthelper_hadronprod helper;
    helper.miny=miny; helper.maxy=maxy;

    helper.N=this; 
    helper.deuteron=deuteron;
    helper.final=final;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.sqrts=sqrts;
    helper.maxpt=maxpt; helper.minpt=minpt;
    helper.ptint=true;

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_yint;
    fun.params=&helper;
    
    // If miny=maxy, don't integrate over y, only over p_T
    if (std::abs(maxy-miny)<0.001)
		return Inthelperf_hadronprod_yint(maxy, &helper);

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(HADRONPROD_YINTPOINTS);

    int status=0; double abserr, result;
    status=gsl_integration_qag(&fun, miny, maxy,
            0, HADRONPROD_INTACCURACY, HADRONPROD_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "Yint failed at " << LINEINFO <<", result " << result <<", relerr "
            << std::abs(abserr/result) << endl;
    }
    
    gsl_integration_workspace_free(workspace);

    return result;
}

double Inthelperf_hadronprod_yint(double y, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    

    // Don't integrate over pt
    if (par->ptint == false)
        return par->N->dHadronMultiplicity_dyd2pt(y, par->pt, 
            par->sqrts, par->frag, par->pdf, par->deuteron, par->final); 
    
    cout << "# yint y=" << y << endl;
    gsl_function fun;
    fun.function=Inthelperf_hadronprod_ptint;
    par->y = y;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(HADRONPROD_PTINTPOINTS);
    status=gsl_integration_qag(&fun, par->minpt, par->maxpt,
            0, HADRONPROD_INTACCURACY, HADRONPROD_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "yint failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << " y " << y << endl;
    }
    
    gsl_integration_workspace_free(workspace);

    return result*2.0*M_PI; // 2\pi from angular integral

}

double Inthelperf_hadronprod_ptint(double pt, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    return pt*par->N->dHadronMultiplicity_dyd2pt(par->y, pt, par->sqrts, par->frag,
        par->pdf, par->deuteron, par->final);   
}

/*
 * Average hadron yield in rapidity range
 * Integrate HadronMultiplicity_dyd2pt over y region and averge
 */
double AmplitudeLib::AverageHadronMultiplicity(double miny, double maxy, double pt, double sqrts, 
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    Inthelper_hadronprod helper;
    helper.N=this; 
    helper.deuteron=deuteron;
    helper.final=final;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.sqrts=sqrts;
    helper.ptint=false;
    helper.pt=pt;

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_yint;
    fun.params=&helper;
    const int YINT_HADRONAVERAGE = 2;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(YINT_HADRONAVERAGE);
    status=gsl_integration_qag(&fun, miny, maxy,
            0, 0.01, YINT_HADRONAVERAGE,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "yint failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << " pt " << pt << endl;
    }
    
    gsl_integration_workspace_free(workspace);
    
    return result / (maxy-miny);
}

/*
 * Douple parton scattering
 */
struct Inthelper_dps
{
    double sqrts;
    FragmentationFunction* fragfun;
    bool deuteron;
    Hadron final;
    AmplitudeLib* N;
    double y1,y2,pt1,pt2;
    double xf1,xf2;
    double z1;
    PDF *pdf;
};


const int DPS_ZINTPOINTS = 6;

double Inthelperf_dps_z1(double z1, void* p);
double Inthelperf_dps_z2(double z2, void* p);

double AmplitudeLib::DPS(double y1, double y2, double pt1, double pt2, double sqrts,
              FragmentationFunction* fragfun, PDF *pdf, bool deuteron, Hadron final)
{	const double ncoll = 15.1;	// central d+Au
	const double C_p = 3;
	Inthelper_dps par;
	par.xf1 = pt1/sqrts*std::exp(y1);
	par.xf2 = pt2/sqrts*std::exp(y2);
	par.y1=y1; par.y2=y2; par.pt1=pt1; par.pt2=pt2; par.N=this;
	par.pdf= pdf; par.sqrts=sqrts; par.deuteron=deuteron; par.final=final;
	par.fragfun=fragfun;
	
	gsl_function fun;
    fun.function=Inthelperf_dps_z1;
    fun.params=&par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_ZINTPOINTS);
    status=gsl_integration_qag(&fun, par.xf1, 1.0,
            0, 0.1, DPS_ZINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "z1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    
    gsl_integration_workspace_free(workspace);
	
	
	result *= (1.0 + 1.0/ncoll) * C_p / std::pow(2.0*M_PI, 4);
	return result;
	
}

double Inthelperf_dps_z1(double z1, void* p)
{
	Inthelper_dps* par = (Inthelper_dps*)p;
	par->z1=z1;
	gsl_function fun;
    fun.function=Inthelperf_dps_z2;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_ZINTPOINTS);
    status=gsl_integration_qag(&fun, par->xf2, 1.0,
            0, 0.1, DPS_ZINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "z2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }

    gsl_integration_workspace_free(workspace);
    return result/SQR(z1);
    
}

double Inthelperf_dps_z2(double z2, void* p)
{
	Inthelper_dps* par = (Inthelper_dps*)p;
	double xp1 = par->xf1 / par->z1;
	double xp2 = par->xf2 / z2;
	double xa1 = xp1 * std::exp(-2.0*par->y1);
	double xa2 = xp2 * std::exp(-2.0*par->y2);
	double ya1 = std::log( par->N->X0() / xa1);
	double ya2 = std::log( par->N->X0() / xa2);
	if (ya1<0 or ya2<0)
	{
		cerr <<"Too large x_A: xa1=" << xa1 <<", xa2=" << xa2 <<", xp1=" << xp1 <<", xp2=" 
			<< xp2 <<", z1=" << par->z1 <<", z2=" << z2 << ", xf1=" << par->xf1
			<<", xf2=" << par->xf2 << " " << LINEINFO <<endl;
		exit(1);
	}
	
	if (xp2 / (1.0-xp1) >= 1.0)
		return 0;
	
	par->N->SetOutOfRangeErrors(false);
	par->N->InitializeInterpolation(ya1);
	double nf1 = par->N->S_k(par->pt1/par->z1, ya1);
	par->N->InitializeInterpolation(ya2);
	double nf2 = par->N->S_k(par->pt2/z2, ya2);
	
	double scale = std::max(par->pt1, par->pt2);
	
	Parton partons[2] = {U, D};
	double xf_d=0;	// x1*x2*f(x1,x2)*D(z1)*D(z2)
	for (int p1ind=0; p1ind<=1; p1ind++)
	{
		for (int p2ind=0; p2ind<=1; p2ind++)
		{
			// Average over proton and neutron
			xf_d += 0.5 * (
				par->pdf->xq(xp1, scale , partons[p1ind])/(xp1)
				* par->pdf->xq( xp2/(1.0-xp1), scale, partons[p2ind] ) / ( xp2/(1.0-xp1) ) 
			   + par->pdf->xq(xp2, scale , partons[p2ind])/(xp2)
				* par->pdf->xq( xp1/(1.0-xp2), scale, partons[p1ind] ) / ( xp1/(1.0-xp2) )
			) * xp1*xp2	
			*  par->fragfun->Evaluate(partons[p1ind], par->final, par->z1, scale)
			 * par->fragfun->Evaluate(partons[p2ind], par->final, z2, scale);	
		}
	}		
	double result = nf1*nf2*xf_d;		
	
	return result/SQR(z2);	
}

/*
 * DPS integrated over pt and y range
 */
struct Inthelper_dpsint
{
	double miny,maxy,minpt,maxpt,sqrts;
	double pt1,y1,y2;
	AmplitudeLib* N;
	FragmentationFunction* fragfun;
	PDF* pdf;
	bool deuteron;
	Hadron final;
};
const int DPS_YINTPOINTS = 1;
const int DPS_PTINTPOINTS = 3;
double Inthelperf_dpsint_y1(double y1, void* p);
double Inthelperf_dpsint_y2(double y2, void* p);
double Inthelperf_dpsint_pt1(double pt1, void* p);
double Inthelperf_dpsint_pt2(double pt2, void* p);
double AmplitudeLib::DPSMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
			FragmentationFunction* fragfun, PDF* pdf, bool deuteron, Hadron final)
{
	cout <<"# DPS contribution integrating yrange " << miny << " - " << maxy 
		<<", ptrange " << minpt << " " << maxpt << endl;
	Inthelper_dpsint par;
	par.N=this;
	par.pdf= pdf; par.sqrts=sqrts; par.deuteron=deuteron; par.final=final;
	par.fragfun=fragfun; par.miny=miny; par.maxy=maxy; par.minpt=minpt;
	par.maxpt=maxpt;
	
	gsl_function fun;
    fun.function=Inthelperf_dpsint_y1;
    fun.params=&par;

	if (std::abs(miny-maxy)<0.1) 
	{
		cout <<"# Not integrating over rapidity range... setting y=" << miny << endl;
		par.y1=miny; par.y2=miny;
		return SQR(2.0*M_PI)*Inthelperf_dpsint_y2(miny, &par); //(2\pi)^2 from angural integrals
	}
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_YINTPOINTS);
    status=gsl_integration_qag(&fun, miny, maxy,
            0, 0.1, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "y1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result*SQR(2.0*M_PI);
	
}
double Inthelperf_dpsint_y1(double y1, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	par->y1=y1;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_y2;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_YINTPOINTS);
    status=gsl_integration_qag(&fun, par->miny, par->maxy,
            0, 0.1, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "y2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}
double Inthelperf_dpsint_y2(double y2, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	cout <<"# y1 " << par->y1 << " y2 "  << y2 << endl;
	par->y2=y2;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_pt1;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);
    status=gsl_integration_qag(&fun, 2, 5,
            0, 0.1, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "pt1int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}

double Inthelperf_dpsint_pt1(double pt1, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	
	par->pt1=pt1;
	gsl_function fun;
	fun.function=Inthelperf_dpsint_pt2;
    fun.params=par;
    
    int status=0; double abserr, result;
    gsl_integration_workspace *workspace 
        = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);
    status=gsl_integration_qag(&fun, 1, pt1,
            0, 0.1, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);

    if (status)
    {
        cerr << "pt2int failed at " << LINEINFO <<", result " << result
            << " relerr " << std::abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(workspace);
	return result;
}

double Inthelperf_dpsint_pt2(double pt2, void* p)
{
	Inthelper_dpsint* par = (Inthelper_dpsint*) p;
	cout <<"# pt1 " << par->pt1 << " pt2 " << pt2 << endl;
	
	return pt2 * par->pt1 * par->N->DPS(par->y1, par->y2, par->pt1, pt2, 
		par->sqrts, par->fragfun, par->pdf, par->deuteron, par->final);
}
