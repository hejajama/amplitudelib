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
    PDF* pdf;
    FragmentationFunction* frag;
};

REAL Inthelperf_hadronprod(REAL z, void *p)
{
    Inthelper_hadronprod* par = (Inthelper_hadronprod*)p;
    REAL x1 = par->xf/z; REAL x2 = x1*std::exp(-2.0*par->y);
    REAL y_A = std::log(0.01/x2);
    par->N->InitializeInterpolation(y_A);
    if (y_A<0)
    {
        cout << "z " << z << " xf " << par->xf << " x1 " << x1 << " x2 " << x2 << " y " << y_A << endl;
        exit(1);
    }

    REAL result = 0;

    // Quark from proton:
    
    // UGD in fundamental representation
    REAL nf = par->N->S_k(par->pt/z, y_A);
    // PDF and fragmentation
    REAL xqf = par->pdf->xq(x1, par->pt, UVAL)*par->frag->Evaluate(U, PI0, z, par->pt)
        + par->pdf->xq(x1, par->pt, DVAL)*par->frag->Evaluate(D, PI0, z, par->pt);

    result = nf*xqf;

    // Adjoint representation, gluon scatters
    nf = par->N->S_k(par->pt/z, y_A, true);
    xqf = par->pdf->xq(x1, par->pt, G)*par->frag->Evaluate(G, PI0, z, par->pt);
    result += nf*xqf;

    return result/SQR(z);
}

REAL AmplitudeLib::dHadronMultiplicity_dyd2pt(REAL y, REAL pt, REAL sqrts)
{
    const REAL K = 1.0; // normalization factor
    // We assume light hadrons
    REAL xf = pt/sqrts*std::exp(y);

    Inthelper_hadronprod helper;
    helper.N=this; helper.y=y; helper.pt=pt; helper.xf=xf;

    KKP fragmentation;
    CTEQ pdf;
    pdf.Initialize();

    helper.pdf=&pdf; helper.frag=&fragmentation;

    REAL result=0; REAL abserr=0;
    const int MULTIPLICITYXINTPOINTS=50;

    
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

// x1: fraction of parent hadron momentum carried by quark,
                        // x2: Bjorken x for the dense system
