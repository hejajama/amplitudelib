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
    REAL y, pt, xf;
    REAL sqrts;
    REAL miny,maxy;
    REAL minpt,maxpt;
    PDF* pdf;
    FragmentationFunction* frag;
    bool deuteron;
    Hadron final;
};

REAL Inthelperf_hadronprod(REAL z, void *p)
{
    // x1: fraction of parent hadron momentum carried by quark,
    // x2: Bjorken x for the dense system
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    REAL x1 = par->xf/z; REAL x2 = x1*std::exp(-2.0*par->y);
    if (x1>1) return 0;
    
    REAL y_A = std::log(par->N->X0()/x2);

    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", z " << z << " xf " <<
            par->xf << " x1 " << x1 << " x2 " << x2 << " y " << y_A << endl;
        return 0;
    }

    bool deuteron = par->deuteron;
    par->N->InitializeInterpolation(y_A);

    REAL result = 0;

    // Quark from proton:
    // UGD in fundamental representation
    REAL nf = par->N->S_k(par->pt/z, y_A);
    // PDF and fragmentation
    REAL xqf = par->pdf->xq(x1, par->pt, U)*par->frag->Evaluate(U, par->final, z, par->pt)
        + par->pdf->xq(x1, par->pt, D)*par->frag->Evaluate(D, par->final, z, par->pt);

    if (deuteron)
    {
        // isospin symmetry, u in p -> d in n
        xqf += par->pdf->xq(x1, par->pt, U)*par->frag->Evaluate(D, par->final, z, par->pt)
        + par->pdf->xq(x1, par->pt, D)*par->frag->Evaluate(U, par->final, z, par->pt);
    }
        
    result = nf*xqf;

    // Adjoint representation, gluon scatters
    REAL na = par->N->S_k(par->pt/z, y_A, true);
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
REAL AmplitudeLib::dHadronMultiplicity_dyd2pt(REAL y, REAL pt, REAL sqrts,
    FragmentationFunction* fragfun, PDF* pdf, bool deuteron, Hadron final )
{
    const REAL K = 1.0; // normalization factor
    // We assume light hadrons
    REAL xf = pt/sqrts*std::exp(y);

    Inthelper_hadronprod helper;
    helper.N=this; helper.y=y; helper.pt=pt; helper.xf=xf;
    helper.deuteron=deuteron;
    helper.final=final;
    helper.pdf=pdf; helper.frag=fragfun;

    REAL result=0; REAL abserr=0;
    const int MULTIPLICITYXINTPOINTS=10;

    
    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_hadronprod;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MULTIPLICITYXINTPOINTS);

    int status;
    status=gsl_integration_qag(&fun, xf, 1.0,
            0, 0.01, MULTIPLICITYXINTPOINTS,
            GSL_INTEG_GAUSS41, workspace, &result, &abserr);
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
REAL Inthelperf_hadronprod_yint(REAL y, void *p);
REAL Inthelperf_hadronprod_ptint(REAL pt, void *p);
const int HADRONPROD_YINTPOINTS=1;
const int HADRONPROD_PTINTPOINTS=2;
const REAL HADRONPROD_INTACCURACY=0.01;

REAL AmplitudeLib::HadronMultiplicity(REAL miny, REAL maxy, REAL minpt, REAL maxpt, REAL sqrts,
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

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_yint;
    fun.params=&helper;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(HADRONPROD_YINTPOINTS);

    int status=0; REAL abserr, result;
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

REAL Inthelperf_hadronprod_yint(REAL y, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    cout << "# yint y=" << y << endl;

    gsl_function fun;
    fun.function=Inthelperf_hadronprod_ptint;
    par->y = y;
    fun.params=par;
    
    int status=0; REAL abserr, result;
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

REAL Inthelperf_hadronprod_ptint(REAL pt, void* p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    return pt*par->N->dHadronMultiplicity_dyd2pt(par->y, pt, par->sqrts, par->frag,
        par->pdf, par->deuteron, par->final);   
}

/*
 * Douple parton scattering
 * dHadronMultiplicity_dyd2pt integrated over miny<y1,y2<maxy
 * pt1>minpt1,  pt1>pt2>mintp2
 * See 1005.4065
 */
struct Inthelper_dps
{
    REAL miny, maxy;
    REAL minpt1, minpt2;
    REAL sqrts;
    FragmentationFunction* fragfun;
    bool deuteron;
    Hadron final;
    AmplitudeLib* N;
    REAL y1,y2,pt1;
    PDF *pdf;
};
REAL Inthelperf_dps_y1(REAL y1, void *p);
REAL Inthelperf_dps_y2(REAL y1, void *p);
REAL Inthelperf_dps_pt1(REAL y1, void *p);
REAL Inthelperf_dps_pt2(REAL y1, void *p);


const int DPS_YINTPOINTS=1;
const int DPS_PTINTPOINTS=2;
const REAL DPS_YINTACCURACY=0.05;
const REAL DPS_PTINTACCURACY=0.05;
const REAL DPS_MAXPT=3.5;

REAL AmplitudeLib::DPS(REAL miny, REAL maxy, REAL minpt1, REAL minpt2, REAL sqrts,
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

    int status=0; REAL abserr, result;
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

REAL Inthelperf_dps_y1(REAL y1, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->y1=y1;
    gsl_function fun;
    fun.function=Inthelperf_dps_y2;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_YINTPOINTS);

    int status=0; REAL abserr, result;
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
REAL Inthelperf_dps_y2(REAL y2, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->y2=y2;
    gsl_function fun;
    fun.function=Inthelperf_dps_pt1;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);

    int status; REAL abserr, result;
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

REAL Inthelperf_dps_pt1(REAL pt1, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;
    par->pt1=pt1;
    gsl_function fun;
    fun.function=Inthelperf_dps_pt2;
    fun.params=par;
    
    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(DPS_PTINTPOINTS);

    int status; REAL abserr, result;
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

REAL Inthelperf_dps_pt2(REAL pt2, void* p)
{
    Inthelper_dps* par = (Inthelper_dps*)p;

    REAL n1,n2;

    //cout << "Evaluating pt1=" << par->pt1 << " pt2=" << pt2 << " y1=" << par->y1
    //<< " y2=" << par->y2 << endl;

            n1 =par->N->dHadronMultiplicity_dyd2pt(par->y1,par->pt1, par->sqrts,
                par->fragfun, par->pdf, par->deuteron, par->final);

            n2 = par->N->dHadronMultiplicity_dyd2pt(par->y2,pt2, par->sqrts,
                par->fragfun, par->pdf, par->deuteron, par->final);
    
    
    return n1*n2*par->pt1*pt2;
}
