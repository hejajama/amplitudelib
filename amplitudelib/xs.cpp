/*
 * AmplitudeLib cross section calculation methods
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
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
    const int MULTIPLICITYXINTPOINTS=10;

    
    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_hadronprod;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MULTIPLICITYXINTPOINTS);

    int status;
    status=gsl_integration_qag(&fun, xf, 1.0,
            0, 0.01, MULTIPLICITYXINTPOINTS,
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
 * dHadronMultiplicity_dyd2pt integrated over miny<y1,y2<maxy
 * pt1>minpt1,  pt1>pt2>mintp2
 * See 1005.4065
 */
struct Inthelper_dps
{
    double miny, maxy;
    double minpt1, minpt2;
    double sqrts;
    FragmentationFunction* fragfun;
    bool deuteron;
    Hadron final;
    AmplitudeLib* N;
    double y1,y2,pt1;
    PDF *pdf;
};
double Inthelperf_dps_y1(double y1, void *p);
double Inthelperf_dps_y2(double y1, void *p);
double Inthelperf_dps_pt1(double y1, void *p);
double Inthelperf_dps_pt2(double y1, void *p);


const int DPS_YINTPOINTS=1;
const int DPS_PTINTPOINTS=2;
const double DPS_YINTACCURACY=0.05;
const double DPS_PTINTACCURACY=0.05;
const double DPS_MAXPT=3.5;

double AmplitudeLib::DPS(double miny, double maxy, double minpt1, double minpt2, double sqrts,
            FragmentationFunction* fragfun, bool deuteron, Hadron final)
{
    Inthelper_dps helper;
    helper.miny=miny; helper.maxy=maxy; helper.minpt1=minpt1; helper.minpt2=minpt2;
    helper.fragfun=fragfun; helper.deuteron=deuteron; helper.final=final;
    helper.sqrts=sqrts;
    helper.N=this;

    CTEQ pdf;
    pdf.Initialize(); helper.pdf=&pdf;
    
    gsl_function fun;
    fun.function=Inthelperf_dps_y1;
    fun.params=&helper;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_YINTPOINTS);

    int status=0; double abserr, result;
    /*status=gsl_integration_qag(&fun, miny, maxy,
            0, DPS_YINTACCURACY, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    */
    gsl_integration_workspace_free(workspace);
    //result = 0.5*(maxy-miny)
    //    * (Inthelperf_dps_y1(miny, &helper) + Inthelperf_dps_y1(maxy, &helper));
    result = (maxy-miny)*Inthelperf_dps_y1(0.5*(maxy+miny), &helper);

    if (status)
        cerr << "y1 integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << endl;
    
    return SQR(2.0*M_PI)*result;    // (2pi)^2 from angular integrals
    
}

double Inthelperf_dps_y1(double y1, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->y1=y1;
    gsl_function fun;
    fun.function=Inthelperf_dps_y2;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_YINTPOINTS);

    int status=0; double abserr, result;
    /*status=gsl_integration_qag(&fun, par->miny, par->maxy,
            0, DPS_YINTACCURACY, DPS_YINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
            */
    gsl_integration_workspace_free(workspace);

    //result = 0.5*(par->maxy-par->miny)
    //    * (Inthelperf_dps_y2(par->miny, par) + Inthelperf_dps_y2(par->maxy, par));
    result = (par->maxy-par->miny) * Inthelperf_dps_y2(0.5*(par->maxy + par->miny), par);

    if (status)
        cerr << "y2 integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << endl;
    cout << "y2 int done" << endl;
    return result;
}
int ptint=0;
double Inthelperf_dps_y2(double y2, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->y2=y2;
    gsl_function fun;
    fun.function=Inthelperf_dps_pt1;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);

    int status; double abserr, result;
    status=gsl_integration_qag(&fun, par->minpt1, DPS_MAXPT,
            0, DPS_PTINTACCURACY, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "pt1 integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << endl;
    cout << "pt1 int done" << endl;
    return result;
}

double Inthelperf_dps_pt1(double pt1, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->pt1=pt1;
    gsl_function fun;
    fun.function=Inthelperf_dps_pt2;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);

    int status; double abserr, result;
    status=gsl_integration_qag(&fun, par->minpt2, pt1,
            0, DPS_PTINTACCURACY, DPS_PTINTPOINTS,
            GSL_INTEG_GAUSS15, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "pt2 integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << endl;
    ptint++;
    cout << "#pt2 int " << ptint << " / " << 15*DPS_PTINTPOINTS << " done" << endl;
    return result*SQR(2.0*M_PI);    //(2\pi)^2: angular part
}

double Inthelperf_dps_pt2(double pt2, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;

    double n1,n2;

    //cout << "Evaluating pt1=" << par->pt1 << " pt2=" << pt2 << " y1=" << par->y1
    //<< " y2=" << par->y2 << endl;

            n1 =par->N->dHadronMultiplicity_dyd2pt(par->y1,par->pt1, par->sqrts,
                par->fragfun, par->pdf, par->deuteron, par->final);

            n2 = par->N->dHadronMultiplicity_dyd2pt(par->y2,pt2, par->sqrts,
                par->fragfun, par->pdf, par->deuteron, par->final);
    
    
    return n1*n2*par->pt1*pt2;
}
