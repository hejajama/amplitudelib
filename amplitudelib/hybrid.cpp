/*
 * Hybrid formalism
 * Part of AmplitudeLib
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include "amplitudelib.hpp"
#include "../tools/tools.hpp"
#include "../pdf/pdf.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "single_inclusive.hpp"
#include <gsl/gsl_integration.h>
#include <cmath>

using namespace Amplitude;

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
    SingleInclusive* xs;
    bool deuteron;
    Hadron final;
    bool ptint;     // if false, don't integrate over pt
    double scale;
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
	
	double scale=par->scale;
	if (scale<0) scale = par->pt;

	
    double y_A = std::log(par->N->X0()/x2);


    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", z " << z << " parton_pt " << par->pt/z << " xf " <<
            par->xf << " x1 " << x1 << " x2 " << x2 << " y_A " << y_A 
            << " y " << par->y << " sqrts " 
            << par->sqrts << " pt " << par->pt << endl ;
        return 0;
        y_A=0;
    }

    bool deuteron = par->deuteron;
    par->N->InitializeInterpolation(x2);
    
    // Quark from proton:
    // UGD in fundamental representation
    double nf = par->N->S_k(par->pt/z, x2);
    // Adjoint representation, gluon scatters
    double na = par->N->S_k(par->pt/z, x2, ADJOINT);

    double result=0;

    // Partons
    for (unsigned int i=0; i<par->xs->Partons().size(); i++)
    {
        if (par->xs->Partons()[i]==LIGHT)
        {
            result =  nf *
                (par->pdf->xq(x1, scale, U)*par->frag->Evaluate(U, par->final, z, scale)
                + par->pdf->xq(x1, scale, D)*par->frag->Evaluate(D, par->final, z, scale)
                + par->pdf->xq(x1, scale, S)*par->frag->Evaluate(S, par->final, z, scale)
                );
            if (deuteron)
            {
                    result += nf * (
                        par->pdf->xq(x1, scale, U)*par->frag->Evaluate(D, par->final, z, scale)
                        + par->pdf->xq(x1, scale, D)*par->frag->Evaluate(U, par->final, z, scale)
                        + par->pdf->xq(x1, scale, S)*par->frag->Evaluate(S, par->final, z, scale)
                        );
            }
        }
        else if (par->xs->Partons()[i]==C or par->xs->Partons()[i]==B)
        {
            double contrib = nf * par->pdf->xq(x1, scale, par->xs->Partons()[i])*par->frag->Evaluate(par->xs->Partons()[i], par->final, z, scale);
            //cout << "C contrib: " << contrib << " pdf " <<  par->pdf->xq(x1, scale, par->xs->Partons()[i]) << " frag " << par->frag->Evaluate(par->xs->Partons()[i], par->final, z, scale) << endl;
            if (deuteron)
                contrib *= 2.0;
            result += contrib;
        }

        else if (par->xs->Partons()[i]==G)
        {
            double contrib = na*par->pdf->xq(x1, scale, G)*par->frag->Evaluate(G, par->final, z, scale);
            if (deuteron)
                contrib *= 2.0;
            result += contrib;
        }
        else
        {
            cerr << "Unknown parton " << par->xs->Partons()[i] << " at " << LINEINFO << endl;
            exit(1);
        }
    }

    return result/SQR(z);
}

/*
 * Calculate dN / dyd^2 pt
 * if deuteron is true (default: false), the probe is deuteron, not proton
 * => use isospin symmetry
 * Default final state particle is PI0
 * If optional parameter scale is given (default: -1), use it when evaluating PDF/FF
 * (used e.g. when calculating DPS contribution (b))
 */
double SingleInclusive::dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
    FragmentationFunction* fragfun, PDF* pdf, Hadron final, bool deuteron, double scale )
{
    // We assume light hadrons
    double xf = pt/sqrts*std::exp(y);
    
    N->SetOutOfRangeErrors(false);
    
    if (xf > 1 or sqrts < 10)
    {
        cerr << "Parameters don't make sense, xf=" << xf << ", sqrts="
            << sqrts << ", y=" << y << " pt=" << pt << " " << LINEINFO << endl;
        return 0;
    }
    Inthelper_hadronprod helper;
    helper.y=y; helper.pt=pt; helper.xf=xf;
    helper.deuteron=deuteron;
    helper.final=final;
    helper.xs = this;
    helper.sqrts=sqrts;
    helper.pdf=pdf; helper.frag=fragfun;
    helper.y=y; helper.scale=scale;
    helper.N=N;

    double result=0; double abserr=0;
    const int MULTIPLICITYXINTPOINTS=40;

    
    gsl_function fun;
    fun.params = &helper;
    fun.function = Inthelperf_hadronprod;

    gsl_integration_workspace *workspace 
     = gsl_integration_workspace_alloc(MULTIPLICITYXINTPOINTS);

    int status;
    status=gsl_integration_qag(&fun, std::max(xf, 0.05), 1.0,
            0, 0.005, MULTIPLICITYXINTPOINTS,
            GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    gsl_integration_workspace_free(workspace);

    if (status)
        cerr << "z integral failed at " << LINEINFO <<": result " << result
        << " relerror " << std::abs(abserr/result) << " y " << y << " pt " << pt << endl;


    result *= 1.0 / SQR(2.0*M_PI);

    return result;
}

/*
 * Set partons included in single inclusive calculations
 */
void SingleInclusive::SetPartons(std::vector<Parton> p)
{
    partons.clear();
    for (unsigned int i=0; i<p.size(); i++)
    {
        partons.push_back(p[i]);
    }
}

std::vector<Amplitude::Parton> & SingleInclusive::Partons()
{
    return partons;
}
