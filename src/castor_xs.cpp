/*
 * Single inclusive cross section calculations using the
 * AmplitudeLib, for Castor kinematics, so total energy is defined
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014-2018
 */

#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/config.hpp"
#include "../tools/tools.hpp"
#include "../pdf/cteq.hpp"
#include "../fragmentation/dss.hpp"
#include "../fragmentation/kkp.hpp"
#include "../fragmentation/hkns.hpp"
#include "../fragmentation/pkhff.hpp"
#include "../pdf/eps09.hpp"
#include "../pdf/ugdpdf.hpp"
#include "../amplitudelib/single_inclusive.hpp"
#include "../amplitudelib/gitsha1.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_integration.h>
#include <ctime>
#include <unistd.h>

using namespace Amplitude;
using namespace std;

const double default_particle_mass = 0.2;
const double castor_min_pseudorapidity = 5.2;
const double castor_max_pseudorapidity = 6.6;
const double rapidity_shift = 0.465;

const int INTWORKSPACEDIVISIONS = 5;
const double INTACCURACY = 0.001;

// Kinematics
// Rapidity from pseudorapidity and pt
double Rapidity(double eta, double pt, double m=default_particle_mass)
{
    return std::log(
                    ( std::sqrt(m*m + std::pow(pt*std::cosh(eta),2)) + pt*std::sinh(eta))
                    / std::sqrt(m*m + pt*pt)
            );
}

// Pseudorapidity from rapidity
double Pseudorapidity(double y, double pt, double m=default_particle_mass)
{
    return std::acosh( std::exp(-y)*std::sqrt( std::pow(std::exp(2.0*y)-1.0,2)*std::pow(m,4) + 2.0*(1+std::exp(4.0*y))*std::pow(m*pt,2) + std::pow(1 + std::exp(2.0*y), 2) * std::pow(pt,4) )
                      / (2.0*pt*std::sqrt( m*m + pt*pt))
                      );
}

double JetEnergy(double y, double pt, double m=default_particle_mass)
{
    return 1.0/2.0 * (std::exp(-y) + std::exp(y)) * std::sqrt(m*m + pt*pt);
}

struct inthelper_castor
{
    SingleInclusive *xs;
    double minE;    // energy bin
    double maxE;
    PDF* pdf;
    double sqrts;
    double y;   // rapidity integrated over
    gsl_integration_workspace* workspace;
    double m;
};

double inthelperf_pt(double pt, void* p);
double inthelperf_y(double y, void* p);

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "# Single parton yield   (c) Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2018 " << endl;
    infostr << "# Latest git commit " << g_GIT_SHA1 << endl;
    infostr << "# Command: ";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);


    if ( argc==1 or (string(argv[1])=="-help" or string(argv[1])=="--help")  )
    {
        cout << "-data [bksolutionfile]" << endl;
        cout << "-minE, -maxE: calorimeter energy range" << endl;
        cout << "-sqrts center-of-mass energy [GeV]" << endl;
        cout << "-x0 val: set x0 value for the BK solutions (overrides the value in BK file)" << endl;
        cout << "-gsl_ft: use GSL to directly calculate Fourier transform" << endl;
        return 0;
        
    }


    double x0=-1;
    double sqrts=5020;
    Order order=LO;
    PDF* pdf=NULL;
    FragmentationFunction* fragfun=NULL;
    std::string datafile="";
    FT_Method ft_method = ACC_SERIES;
    std::string datafile_probe="";
    bool deuteron=false;
    double ptstep=0.1;
    double minE=-1; double maxE=-1;

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-x0")
            x0 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile=argv[i+1];
        else if (string(argv[i])=="-sqrts")
            sqrts = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-minE")
            minE  = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxE [default: minE + 100]")
            maxE  = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-lo")
            order=LO;
        else if (string(argv[i])=="-nlo")
            order=NLO;
        else if (string(argv[i])=="-pdf")
        {
			if (string(argv[i+1])=="cteq")
				pdf = new CTEQ();
			else if (string(argv[i+1])=="eps09")
            {
				pdf = new EPS09();
				pdf->SetA(StrToInt(argv[i+2]));
			}
            else if (string(argv[i+1])=="ugd")
			{
				AmplitudeLib* ugdn = new AmplitudeLib(string(argv[i+2]));
				pdf = new UGDPDF(ugdn, StrToReal(argv[i+3]));
			}
			else
			{
				cerr << "Unknown PDF type " << argv[i+1] << endl;
				return -1;
			}
		}
        else if (string(argv[i])=="-gsl_ft")
            ft_method = GSL;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    if (maxE < 0)
        maxE  = minE + 100;
    // Default PDF and FF
    if (pdf==NULL)
        pdf = new CTEQ();
    
    // Read data
    AmplitudeLib N(datafile);

    N.SetFTMethod(ft_method);
    SingleInclusive xs(&N);

    
    
    
    if (x0>0)
    {
        N.SetX0(x0);
    }
   

    pdf->SetOrder(order);
    

    time_t now = time(0);
    string today = ctime(&now);
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# AmplitudeLib v. " << N.Version()  << " running on " << hostname << endl;
    cout <<"# Now is " << today ;
    cout <<"#"<<endl;
	delete[] hostname;

    cout << "# " << N.GetString() << endl;
    cout << "# PDF: " << pdf->GetString() << endl;
    // Print quarks
    std::vector<Parton> ps;
    ps.push_back(U);
	ps.push_back(UBAR);
	ps.push_back(D);
	ps.push_back(DBAR);
	ps.push_back(S);
	ps.push_back(SBAR);
    ps.push_back(C);
	ps.push_back(CBAR);
	
    ps.push_back(G);
    xs.SetPartons(ps);
    cout <<"# Partons: " ;
    for (unsigned int i=0; i<xs.Partons().size(); i++)
    {
        cout << PartonToString(xs.Partons()[i]) << " ";
    }
    cout << endl;

    //////////////////////////////////////
    inthelper_castor par;
    par.xs = &xs;
    par.minE=minE;
    par.maxE=maxE;
    par.sqrts=sqrts;
    par.pdf=pdf;
    par.m=default_particle_mass;
    
    gsl_function fun;
    fun.params=&par;
    fun.function = inthelperf_y;
    
    gsl_integration_workspace *workspace
    = gsl_integration_workspace_alloc(INTWORKSPACEDIVISIONS);
    gsl_integration_workspace *workspace_internal
    = gsl_integration_workspace_alloc(INTWORKSPACEDIVISIONS);
    par.workspace = workspace_internal;
    int status; double result,abserr;
    double miny=castor_min_pseudorapidity-2; double maxy = castor_max_pseudorapidity+2;
    status = gsl_integration_qag(&fun, miny, maxy, 0, INTACCURACY, INTWORKSPACEDIVISIONS, GSL_INTEG_GAUSS51, workspace, &result, &abserr);
    if (status)
        cerr << "y integral failed at " << LINEINFO <<": result " << result
            << " relerror " << std::abs(abserr/result) << endl;
    
    gsl_integration_workspace_free(workspace);
    gsl_integration_workspace_free(workspace_internal);
    
    
    cout << "#" << endl;
    cout << "# minE   maxE   yield [note: for p, multiply by sigma0/2 / sigma_inel, for nuke: do b int]" << endl;
    cout << minE << " " << maxE << " " << result << endl;

    

    return 0;
}

double inthelperf_y(double y, void* p)
{
    inthelper_castor *par = (inthelper_castor*)p;
    par->y=y;
    gsl_function fun;
    fun.params=par;
    fun.function = inthelperf_pt;
    
    // If you use this, then remember to free it
    //gsl_integration_workspace *workspace
    //= gsl_integration_workspace_alloc(INTWORKSPACEDIVISIONS);
    int status; double result,abserr;
    
    status = gsl_integration_qag(&fun, 0.2, 20, 0, INTACCURACY, INTWORKSPACEDIVISIONS, GSL_INTEG_GAUSS51, par->workspace, &result, &abserr);
    
    if (status)
        cerr << "pt integral failed, result " << result << " relerr " << std::abs(abserr/result) << " y " << y << endl;
    
    return result;
}

double inthelperf_pt(double pt, void* p)
{

    
    inthelper_castor *par = (inthelper_castor*)p;
    
    
 
    // Check kinematics
    double energy = JetEnergy(par->y, pt);
    
    if (energy < par->minE or energy > par->maxE)
    {
        //cout << "Out of kinematics " << par->y << " pt " << pt << endl;
        //cout << par->y << " " << pt << endl;
        return 0;
    }
    
    double eta = Pseudorapidity(par->y, pt);
    if (eta+rapidity_shift < castor_min_pseudorapidity or eta+rapidity_shift > castor_max_pseudorapidity)
    {
        //cout << "Out of castor acceptance  " << par->y << " pt " << pt << endl;
        //cout << par->y << " " << pt << endl;
        return 0;
    }
    
    // integration measure: we want to compute
    // \int dEdy xs(pt,y) = \int dpt dy *Jacobian* xs(pt,y)
    // The Jacobian is dE/dpt = (exp(-y)+exp(y))*pt / (2*sqrt(m^2+pt^2))
    double jacobian = (std::exp(-par->y) + std::exp(par->y))*pt / (2.0 * std::sqrt(par->m*par->m + pt*pt));
    
    // last parameters: false=no deuteron, -1 = scale is pt
    return jacobian*par->xs->dHadronMultiplicity_dyd2pt_parton(par->y, pt, par->sqrts,
                                                      par->pdf, false,  -1 );
}

    
