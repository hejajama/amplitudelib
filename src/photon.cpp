#include "photon.hpp"
#include <cmath>
#include <iostream>
#include "../amplitudelib/amplitudelib.hpp"
#include <gsl/gsl_integration.h>
#include "../tools/tools.hpp"
#include <gsl/gsl_errno.h>
#include "../pdf/pdf.hpp"
#include "../pdf/cteq.hpp"

using namespace std;

int INTERVALS_PHOTON_PT = 4;
int INTERVALS_PHOTON_Y = 4;
int INTERVALS_PHOTON_PHI = 4;
const double INTACCURACY_PHOTON = 0.001;

const double MINY = -5;
const double MAXY = 5;
const double MINPT = 0.1;
const double MAXPT = 20;

const double MAX_XP = 0.9;

double PartonCharge(Amplitude::Parton p)
{
    switch (p)
    {
        case U:
        case C:
            return 2.0/3.0;
        case D:
        case S:
            return -1.0/3.0;
        default:
            cerr <<"Unknown parton!" << endl;
            return 0;
        
    }
    return 0;
}

/*
 * Integrated over quark rapidity and momentum and azimuthal angle
 *  alpha_em/(2pi)^3 is outside this
 */
double IsolatedPhoton::DifferentialPhotonCrossSection(double k, double y_k, double l, double y_l, double phi)

{
    // Check isolation
    double  dphi = min(phi, 2.0*M_PI-phi);
    double dsqr = pow(y_k - y_l, 2.0) + pow(dphi, 2.0);
    
    //cout <<
    if (dsqr < isolation_cut*isolation_cut)
        return 0;
    
    // x in probe
    double xp = k/sqrts * exp(y_k) + l/sqrts*exp(y_l);
    double xg = k/sqrts * exp(-y_k) + l/sqrts*exp(-y_l);
    
    // Check if kinematically allowed
    if (xp > MAX_XP or xg > 0.01 )
        return 0;
    
    double z = k / (xp * sqrts) * exp(y_k); // Amir, Jamal, (A15)
    
    double integrand = (k*k + l*l + 2.0*k*l*cos(phi)) / ( pow(z*l, 2.0) + pow((1.0-z)*k, 2.0) - 2.0*z*(1.0-z)*k*l*cos(phi));
    
    double k_plus_l = sqrt(k*k + l*l + 2.0*l*k*cos(phi));
    
    N->InitializeInterpolation(xg);
    
    integrand *= N->S_k(k_plus_l, xg);
    
    if (isnan(integrand) or isinf(integrand))
    {
        //cerr << integrand  << "z: " << z << " l: " << l << " y_l " << y_l << " k " << k << " y_k " << y_k << " phi " << phi << " dist " << sqrt(dsqr) << endl;
    }
    
    integrand *= z*z*(1.0 + pow(1.0-z,2.0))/(k*k);
    
    integrand *= l; //Jacobian
    
    integrand /= M_PI; // (8) is for dsigma/dkt^2, want dsigma/d^2 kt
    
    return integrand;
    
}


struct Inthelper_isolated_photon
{
    IsolatedPhoton* photon;
    CTEQ* pdf;
    double y_k;
    double k;
    double y_l;
    double l;
    double sqrts;
    gsl_integration_workspace *workspace_ltint;
    gsl_integration_workspace *workspace_phiint;
};

double inthelperf_photon_phi(double phi, void* p)
{
    Inthelper_isolated_photon* par = (Inthelper_isolated_photon*) p;
    return par->photon->DifferentialPhotonCrossSection(par->k, par->y_k, par->l, par->y_l, phi);
}

double inthelperf_photon_lt(double l, void* p)
{
    Inthelper_isolated_photon* par = (Inthelper_isolated_photon*) p;
    par->l = l;
    gsl_function fun;
    
    double xp =par->k/par->sqrts * exp(par->y_k) + par->l/par->sqrts*exp(par->y_l);
    if (xp > MAX_XP)
        return 0;
    
    fun.params=par; fun.function = inthelperf_photon_phi;
    double result,abserr;
    int status = gsl_integration_qag(&fun, 0, M_PI, 0, INTACCURACY_PHOTON,
                                     INTERVALS_PHOTON_PHI, GSL_INTEG_GAUSS15, par->workspace_phiint, &result, &abserr);
    
    
    
    if (status)
    {
        //cerr << "phi integral result " << result << " relerr " << abs(abserr/result) << endl;
    }
    
    result *= 2.0; // As we integrate [0,pi] here
    
    // PDF
    double partonsum=0;
    
    
    double scale = max(par->k,par->l);
    vector<Amplitude::Parton> partons = par->photon->GetPartons();
    for (unsigned int i=0; i<partons.size(); i++)
    {
        double eq = abs(PartonCharge(partons[i]));
        partonsum += eq*eq*par->pdf->xq(xp, scale, partons[i])/xp;
    }

    
    return result*partonsum;
}

double inthelperf_photon_rapidity(double y_l, void* p)
{
    Inthelper_isolated_photon* par = (Inthelper_isolated_photon*) p;
    par->y_l = y_l;
    gsl_function fun;
    
    fun.params=par; fun.function = inthelperf_photon_lt;
    double result,abserr;
    int status = gsl_integration_qag(&fun, MINPT, MAXPT, 0, INTACCURACY_PHOTON,
                    INTERVALS_PHOTON_PT, GSL_INTEG_GAUSS15, par->workspace_ltint, &result, &abserr);
    
    if (status)
    {
        //cerr << "pt integral result " << result << " relerr " << abs(abserr/result) << endl;
    }
    
    return result;
}


double IsolatedPhoton::PhotonCrossSection(double k, double y_k)
{
    Inthelper_isolated_photon helper;
    helper.photon=this;
    helper.k = k;
    helper.y_k  = y_k;
    helper.sqrts=sqrts;
    
    CTEQ pdf; pdf.SetOrder(Amplitude::LO);
    helper.pdf = &pdf;
    
    // Initialize integration workspace for l_t integral here, not every time
    helper.workspace_ltint = gsl_integration_workspace_alloc(INTERVALS_PHOTON_PT);
    helper.workspace_phiint = gsl_integration_workspace_alloc(INTERVALS_PHOTON_PHI);
    gsl_integration_workspace *ws = gsl_integration_workspace_alloc(INTERVALS_PHOTON_Y);
    
    gsl_function fun; fun.params=&helper;
    fun.function=inthelperf_photon_rapidity;
    
    double result,abserr;
    
    
    int status = gsl_integration_qag(&fun, MINY, MAXY, 0, INTACCURACY_PHOTON,
                       INTERVALS_PHOTON_Y, GSL_INTEG_GAUSS15, ws, &result, &abserr);
    
    if (status)
    {
        //cerr << "rapidity integral result " << result << " relerr " << abs(abserr/result) << endl;
    }
    gsl_integration_workspace_free(ws);
    gsl_integration_workspace_free(helper.workspace_ltint);
    gsl_integration_workspace_free(helper.workspace_phiint);
    
    // Constants
    result *= ALPHA_e/pow(2.0*M_PI, 3.0);

    return result;
}




int main(int argc, char* argv[])
{
    gsl_set_error_handler(&ErrHandler);
    string datafile;
    double R = 0.1;
    double y = 3;
    bool gsl_ft = false;
    double sqrts=5020;
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-data")
            datafile = argv[i+1];
        else if (string(argv[i])=="-R")
            R = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-sqrts")
            sqrts= StrToReal(argv[i+1]);
        else if (string(argv[i])=="-y")
            y= StrToReal(argv[i+1]);
        else if (string(argv[i])=="-gsl_ft")
            gsl_ft=true;
        else if (string(argv[i]) == "-integration_intervals")
        {
            int intervals = StrToInt(argv[i+1]);
            INTERVALS_PHOTON_PT=intervals;
            INTERVALS_PHOTON_PHI=intervals;
            INTERVALS_PHOTON_Y = intervals;
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
        
    }
    AmplitudeLib N(datafile);
    if (gsl_ft)
        N.SetFTMethod(Amplitude::GSL);
    IsolatedPhoton photon(&N);
    photon.SetIsolationCut(R);
    photon.SetSqrts(sqrts);
    
    cout << "# Datafile " << datafile << " photon rapidity y=" << y << " isolation cut " << R << " sqrts=" << sqrts << " GeV" << endl;
    if (gsl_ft) cout << "# FT method: GSL " << endl;
    else cout << "# FT method: Acc. series" << endl;
    cout << "# Integration intervals: " << INTERVALS_PHOTON_PT << endl;
    
    for (double pt=1; pt<10; pt+=0.5)
    {
        cout << pt << " " << photon.PhotonCrossSection(pt, y) << endl;
    }
    return 0;
}
