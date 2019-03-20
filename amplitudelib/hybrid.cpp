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
		cerr << "Entering kinematically fobidden region at y=" << par->y <<", pt=" << par->pt << ", z=" << z << " " << LINEINFO << endl;
		return 0;
	}
	
	double scale=par->scale;
	if (scale<0) scale = par->pt;

	if (x2 > par->N->X0()) 
	{
			// We can only end up here if y < 2.3
			x2 = par->N->X0();
	}
    double y_A = std::log(par->N->X0()/x2);


    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO <<", z " << z << " parton_pt " << par->pt/z << " xf " <<
            par->xf << " x1 " << x1 << " x2 " << x2 << " y_A " << y_A 
            << " y " << par->y << " sqrts " 
            << par->sqrts << " pt " << par->pt << endl ;
        return 0;
    }
    bool deuteron = par->deuteron;
    par->N->InitializeInterpolation(x2);
    
    // Quark from proton:
    // UGD in fundamental representation
    double nf = par->N->S_k(par->pt/z, x2);
    // Adjoint representation, gluon scatters
    double na = par->N->S_k(par->pt/z, x2, ADJOINT);

	if (nf < 0) nf = 0;
	if (na < 0) na = 0;

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
        else if (par->xs->Partons()[i]!=G)
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
	double minz = std::max(xf, 0.05);
    status=gsl_integration_qag(&fun, minz, 1.0,
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
 * Hadron production in parton level
 */
double SingleInclusive::dHadronMultiplicity_dyd2pt_parton(double y, double pt, double sqrts,
                                                          PDF* pdf, bool deuteron,  double scale )
{
    // x1: fraction of parent hadron momentum carried by quark,
    // x2: Bjorken x for the dense system
    double x1 = pt/sqrts*std::exp(y); double x2 = x1*std::exp(-2.0*y);
    if (x1>1)
    {
        cerr << "Entering kinematically fobidden region at y=" << y <<", pt=" << pt << " " << LINEINFO << endl;
        return 0;
    }
    
    if (scale<0) scale = pt;
    
    if (x2 > N->X0())
        x2 = N->X0();
    
    double y_A = std::log(N->X0()/x2);
    
    
    if (y_A<0)
    {
        cerr << "Negative rapidity at " << LINEINFO << " parton_pt " << pt <<
        " x1 " << x1 << " x2 " << x2 << " y_A " << y_A
        << " y " << y << " sqrts "
        << sqrts << " pt " << pt << endl ;
        return 0;
    }
    N->InitializeInterpolation(x2);
    
    // Quark from proton:
    // UGD in fundamental representation
    double nf = N->S_k(pt, x2);
    // Adjoint representation, gluon scatters
    double na = N->S_k(pt, x2, ADJOINT);
    
    if (nf < 0) nf = 0;
    if (na < 0) na = 0;
    
    double result=0;
    
    // Partons
    for (unsigned int i=0; i<Partons().size(); i++)
    {
        if (Partons()[i]==LIGHT)
        {
            result =  nf *
            (pdf->xq(x1, scale, U)
             + pdf->xq(x1, scale, D)
             + pdf->xq(x1, scale, S)
             );
            if (deuteron)
            {
                result *= 2.0;
            }
        }
        else if (Partons()[i]!=G)
        {
            double contrib = nf * pdf->xq(x1, scale, Partons()[i]);

            if (deuteron)
                contrib *= 2.0;
            result += contrib;
        }
        
        else if (Partons()[i]==G)
        {
            double contrib = na*pdf->xq(x1, scale, G);
            if (deuteron)
                contrib *= 2.0;
            result += contrib;
        }
        else
        {
            cerr << "Unknown parton " << Partons()[i] << " at " << LINEINFO << endl;
            exit(1);
        }
    }
    
    return result/SQR(2.0*M_PI);
    
    
}


/*
 * Parton level double parton scattering
 * Compute 1/(2pi)^4 x_1 x_2 D(x_1,x_2) S(pt1)S(pt2)
 * Summed over all parton species
 */
double SingleInclusive::dHadronMultiplicity_dyd2pt_parton_dps(double y1, double pt1, double y2, double pt2, double sqrts,
                                         PDF* pdf, bool deuteron,  double scale )
{
    
    
    // xp: bjorken x for probe
    // xA: bjorken x for target
    double xp1 = pt1/sqrts*std::exp(y1); double xA1 = xp1*std::exp(-2.0*y1);
    double xp2 = pt2/sqrts*std::exp(y2); double xA2 = xp1*std::exp(-2.0*y2);
    if (xp1>1 or xp2>1 or xp1+xp2>1)
    {
        cerr << "Entering kinematically fobidden region at y1=" << y1 <<", pt1=" << pt1 << " y2=" << y2 << " pt2=" << pt2 << LINEINFO << endl;
        return 0;
    }
    
    if (scale<0) scale = (pt1+pt2)/2.0;
    
    
    
    if (xA1 > N->X0() or xA2 > N->X0() )
    {
        cerr << "Negative rapidity at " << LINEINFO << " xA1 " << xA1 <<
        " xA2 " << xA2 << " sqrts "
        << sqrts << endl ;
        return 0;
    }
   
    
    // Quark from proton:
     N->InitializeInterpolation(xA1);
    // UGD in fundamental representation
    double nf1 = N->S_k(pt1, xA1);
    // Adjoint representation, gluon scatters
    double nA1 = N->S_k(pt1, xA1, ADJOINT);
    

    N->InitializeInterpolation(xA2);
    double nf2 = N->S_k(pt2, xA2);
    double nA2 = N->S_k(pt2, xA2, ADJOINT);
    
    if (nf1 < 0) nf1 = 0;
    if (nA1 < 0) nA1 = 0;
    if (nf2 < 0) nf2 = 0;
    if (nA2 < 0) nA2 = 0;
    
    double scale_1 = scale;
    double scale_2 = scale;
    if (scale < 0)  // Automatically set scale
    {
        scale_1 = pt1; scale_2=pt2;
        if (scale_1 < pdf->MinQ()) scale_1 = pdf->MinQ();
        if (scale_2 < pdf->MinQ()) scale_2 = pdf->MinQ();
    }
    
    double result=0;
    
    if (deuteron)
    {
        cerr << "Deuteron is not supported in DPS calculation! " << LINEINFO << endl;
        exit(1);
    }
    
    // Partons
    for (unsigned int i=0; i<Partons().size(); i++)
    {
        for (unsigned int j=0; j<Partons().size(); j++)
        {
            double sa_1=0;
            double sa_2=0;
            if (Partons()[i]==G)
                sa_1 = nA1;
            else if (Partons()[i]!=G and Partons()[i]!=LIGHT)
                sa_1 = nf1;
            else if (Partons()[i]==LIGHT)
            {
                cerr << "Light partons not supported in DPS" << endl;
                exit(1);
            }
            else
            {
                cerr << "WTF parton " <<Partons()[i] << " at " << LINEINFO << endl;
                exit(1);
            }
            
            if (Partons()[j]==G)
                sa_2 = nA2;
            else if (Partons()[j]!=G and Partons()[j]!=LIGHT)
                sa_2 = nf2;
            else if (Partons()[j]==LIGHT)
            {
                cerr << "Light partons not supported in DPS" << endl;
                exit(1);
            }
            else
            {
                cerr << "WTF parton " <<Partons()[i] << " at " << LINEINFO << endl;
                exit(1);
            }
            
            double f_i = pdf->xq(xp1, scale_1, Partons()[i])/xp1;
            double f_i_scaled = pdf->xq(xp1/(1.0-xp2), scale_1, Partons()[i]) / (xp1/(1.0-xp2));
            double f_j = pdf->xq(xp2, scale_2, Partons()[j])/xp2;
            double f_j_scaled = pdf->xq(xp2/(1.0-xp1), scale_2, Partons()[j]) / (xp2/(1.0-xp1));
            
            double dpdf = 0.5*xp1*xp2*( f_i * f_j_scaled + f_i_scaled*f_j );
            // For testing: this should reproduce the factorized result
            //dpdf = xp1*xp2*f_i*f_j;
            
            result += sa_1*sa_2 *dpdf;
        }
            
            
            
        
    }
    
    return result/std::pow(2.0*M_PI, 4.0);
    
    
}

/*
 * Parton level triple parton scattering
 * Compute 1/(2pi)^6 x_1 x_2 x_3 D(x_1,x_2,x_3) S(pt1)S(pt2)S(pt3)
 * Summed over all parton species
 */
double SingleInclusive::dHadronMultiplicity_dyd2pt_parton_3ps(double y1, double pt1, double y2, double pt2, double y3, double pt3, double sqrts,    PDF* pdf, bool deuteron,  double scale )
{
    
    
    // xp: bjorken x for probe
    // xA: bjorken x for target
    double xp1 = pt1/sqrts*std::exp(y1); double xA1 = xp1*std::exp(-2.0*y1);
    double xp2 = pt2/sqrts*std::exp(y2); double xA2 = xp1*std::exp(-2.0*y2);
    double xp3 = pt3/sqrts*std::exp(y3); double xA3 = xp1*std::exp(-2.0*y3);
    if (xp1>1 or xp2>1 or xp3>1 or xp1+xp2+xp3>1)
    {
        cerr << "Entering kinematically fobidden region at y1=" << y1 <<", pt1=" << pt1 << " y2=" << y2 << " pt2=" << pt2 << LINEINFO << endl;
        return 0;
    }
    
    if (xA1 > N->X0() or xA2 > N->X0() or xA3 > N->X0() )
    {
        cerr << "Negative rapidity at " << LINEINFO << " xA1 " << xA1 <<
        " xA2 " << xA2 << " sqrts "
        << sqrts << endl ;
        return 0;
    }
    
    
    // Quark from proton:
    N->InitializeInterpolation(xA1);
    // UGD in fundamental representation
    double nf1 = N->S_k(pt1, xA1);
    // Adjoint representation, gluon scatters
    double nA1 = N->S_k(pt1, xA1, ADJOINT);
    
    
    N->InitializeInterpolation(xA2);
    double nf2 = N->S_k(pt2, xA2);
    double nA2 = N->S_k(pt2, xA2, ADJOINT);
    
    N->InitializeInterpolation(xA3);
    double nf3 = N->S_k(pt3, xA3);
    double nA3 = N->S_k(pt3, xA3, ADJOINT);
    
    if (nf1 < 0) nf1 = 0;
    if (nA1 < 0) nA1 = 0;
    if (nf2 < 0) nf2 = 0;
    if (nA2 < 0) nA2 = 0;
    if (nf2 < 0) nf3 = 0;
    if (nA2 < 0) nA3 = 0;
    
    double scale_1 = scale;
    double scale_2 = scale;
    double scale_3 = scale;
    if (scale < 0)  // Automatically set scale
    {
        scale_1 = pt1; scale_2=pt2; scale_3 = pt3;
        if (scale_1 < pdf->MinQ()) scale_1 = pdf->MinQ();
        if (scale_2 < pdf->MinQ()) scale_2 = pdf->MinQ();
        if (scale_3 < pdf->MinQ()) scale_3 = pdf->MinQ();
    }
    
    double result=0;
    
    if (deuteron)
    {
        cerr << "Deuteron is not supported in DPS calculation! " << LINEINFO << endl;
        exit(1);
    }
    
    // Partons
    for (unsigned int i=0; i<Partons().size(); i++)
    {
        for (unsigned int j=0; j<Partons().size(); j++)
        {
            for (unsigned int k=0; k<Partons().size(); k++)
            {
                double sa_1=0;
                double sa_2=0;
                double sa_3=0;
                
                if (Partons()[i]==G)
                    sa_1 = nA1;
                else if (Partons()[i]!=G and Partons()[i]!=LIGHT)
                    sa_1 = nf1;
                else if (Partons()[i]==LIGHT)
                {
                    cerr << "Light partons not supported in DPS" << endl;
                    exit(1);
                }
                else
                {
                    cerr << "WTF parton " <<Partons()[i] << " at " << LINEINFO << endl;
                    exit(1);
                }
                
                if (Partons()[j]==G)
                    sa_2 = nA2;
                else if (Partons()[j]!=G and Partons()[j]!=LIGHT)
                    sa_2 = nf2;
                else if (Partons()[j]==LIGHT)
                {
                    cerr << "Light partons not supported in DPS" << endl;
                    exit(1);
                }
                else
                {
                    cerr << "WTF parton " <<Partons()[j] << " at " << LINEINFO << endl;
                    exit(1);
                }
                
                if (Partons()[k]==G)
                    sa_3 = nA3;
                else if (Partons()[k]!=G and Partons()[k]!=LIGHT)
                    sa_3 = nf3;
                else if (Partons()[k]==LIGHT)
                {
                    cerr << "Light partons not supported in DPS" << endl;
                    exit(1);
                }
                else
                {
                    cerr << "WTF parton " <<Partons()[k] << " at " << LINEINFO << endl;
                    exit(1);
                }
                
                /////TODO TÄSTÄ ETEENPÄIN
                double f_i = pdf->xq(xp1, scale_1, Partons()[i])/xp1;
                double f_i_scaled = pdf->xq(xp1/(1.0-xp2 - xp3), scale_1, Partons()[i]) / (xp1/(1.0-xp2-xp3));
                double f_j = pdf->xq(xp2, scale_2, Partons()[j])/xp2;
                double f_j_scaled = pdf->xq(xp2/(1.0-xp1-xp3), scale_2, Partons()[j]) / (xp2/(1.0-xp1-xp3));
                double f_k = pdf->xq(xp3, scale_3, Partons()[j])/xp3;
                double f_k_scaled = pdf->xq(xp3/(1.0-xp1-xp2), scale_3, Partons()[k]) / (xp3/(1.0-xp1-xp2));
                
                double dpdf = 1.0/3.0*xp1*xp2*xp3*( f_i * f_j * f_k_scaled + f_i * f_j_scaled * f_k
                                                        + f_i_scaled*f_j*f_k );
                // For testing: this should reproduce the factorized result
                //dpdf = xp1*xp2*f_i*f_j;
                
                result += sa_1*sa_2*sa_3 *dpdf;
            }
        
        }
        
        
    }
    
    return result/std::pow(2.0*M_PI, 6.0);
    
    
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
