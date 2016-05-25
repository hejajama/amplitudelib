/*
 * Single inclusive cross section calculations using the
 * AmplitudeLib
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
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
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <unistd.h>

using namespace Amplitude;
using namespace std;

enum Mode
{
    HYBRID_PT,
    KTFACT_PT,
    HYBRID_PT_AVG,
    HYBRID_PARTON,
    KTFACT_PARTON,
    HYBRID_MULTIPLICITY
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "# Single inclusive yield calculator   (c) Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014 " << endl;
    infostr << "# Command: ";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);


    if ( argc==1 or (string(argv[1])=="-help" or string(argv[1])=="--help")  )
    {
		cout << "==== USAGE ====" << endl;
        cout << "====== Hybrid formalism ======" <<endl;
        cout << "-pt_spectrum p1 p2: compute differential particle p2 (pi0, ch) production yield using hybrid formalism, probe is proton (p) or deuteron (d)" << endl;
        cout << "-pt_spectrum_avg: same as above, but average over y region, must set miny and maxy" << endl;
        cout << "-hadronprod p1 p2: integrated over pt and y range (defined using -miny, -maxy, -minpt, -maxpt)" << endl;
        cout << "-pt_spectrum_parton: parton devel differential yield using the hybrid formalism" << endl;
        cout << endl;
        cout << "====== Kt-factorization ======" <<endl;
        cout << "-pt_spectrum_ktfact p: differential particle p production yield using the kt-factorization " << endl;
        cout << "-pt_spectrum_ktfact_parton: differential gluon production yield using the kt-factorization" << endl;
        cout << "-ktfact_probe datafile: in asymmetric collisions dipole amplitude for the probe " << endl;
        cout << "-fixed_alphas: use fixed coupling in kt-factorization" << endl;
        cout << endl;
        cout << "====== General ======" << endl;
        cout << "-data bksolutionfile"<<endl;
        cout << "-y rapidity";
        cout << "-sqrts center-of-mass energy [GeV]" << endl;
        cout << "-minpt, -maxpt, -miny, -maxy" << endl;
        cout << "-ptstep stepsize" << endl;
        cout << "-pdf pdf: possible pdfs: cteq ugd [bkfile sigma0]" << endl;
        cout << "-fragfun ff" << endl;
        cout << "-nlo/-lo: use NLO/LO distributions" << endl;
        cout << "-x0 val: set x0 value for the BK solutions (overrides the value in BK file)" << endl;
        cout << "-gsl_ft: use GSL to directly calculate Fourier transform" << endl;
        return 0;
        
    }


    double x0=-1;
    double y=-1;
    double sqrts=0;
    Order order=LO;
    PDF* pdf=NULL;
    FragmentationFunction* fragfun=NULL;
    std::string datafile="";
    FT_Method ft_method = ACC_SERIES;
    std::string datafile_probe="";
    Mode mode=HYBRID_PT;
    Hadron final_particle=PI0;
    bool deuteron=false;
    double ptstep=0.1;
    double minpt=1, maxpt=2;
    double miny=0; double maxy=1;

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-x0")
            x0 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile=argv[i+1];
        else if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-sqrts")
            sqrts = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxpt")
            maxpt = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-minpt")
            minpt = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-lo")
            order=LO;
        else if (string(argv[i])=="-nlo")
            order=NLO;
        else if (string(argv[i])=="-pt_spectrum" or string(argv[i])=="-pt_spectrum_avg")
        {
            if (string(argv[i])=="-pt_spectrum")
                mode=HYBRID_PT;
            else
                mode=HYBRID_PT_AVG;
            if (string(argv[i+1])=="p")
                deuteron=false;
            else if (string(argv[i+1])=="d")
                deuteron=true;
            else
            {
                cerr << "Invalid probe particle type " << argv[i+1] << endl;
                exit(1);
            }

            if (string(argv[i+2])=="pi0")
                final_particle = PI0;
            else if (string(argv[i+2])=="ch")
                final_particle = H; // charged hadrons
            else if (string(argv[i+2])=="hm")    // negative hadrons
                final_particle = HM;
            else
            {
                cerr << "Invalid final state particle " << argv[i+2] << endl;
                exit(1);
            }
        }
        else if (string(argv[i])=="-pt_spectrum_parton")
			mode=HYBRID_PARTON;
        else if (string(argv[i])=="-pt_spectrum_ktfact")
		{
			mode = KTFACT_PT;
			if (string(argv[i+1])=="pi0")
                final_particle = PI0;
            else if (string(argv[i+1])=="ch")
                final_particle = H; // charged hadrons
            else if (string(argv[i+1])=="hm")    // negative hadrons
                final_particle = HM;
            else
			{
				cerr << "Unknown final particle " << argv[i+1] << endl;
				return -1;
			}
        }
        else if (string(argv[i])=="-pt_spectrum_ktfact_parton")
        {
			mode=KTFACT_PARTON;
		}
        
        else if (string(argv[i])=="-ktfact_probe")
        {
			datafile_probe=argv[i+1];
		}
			
        else if (string(argv[i])=="-hadronprod_int")
        {
            mode=HYBRID_MULTIPLICITY;
            if (string(argv[i+1])=="p")
                deuteron=false;
            else if (string(argv[i+1])=="d")
                deuteron=true;
            else
            {
                cerr << "Invalid probe particle type " << argv[i+1] << endl;
                exit(1);
            }

            if (string(argv[i+2])=="pi0")
                final_particle = PI0;
            else if (string(argv[i+2])=="ch")
                final_particle = H; // charged hadrons
            else if (string(argv[i+2])=="hm")    // negative hadrons
                final_particle = HM;
            else
            {
                cerr << "Invalid final state particle " << argv[i+2] << endl;
                exit(1);
            }
        }
        else if (string(argv[i])=="-ptstep")
            ptstep = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-fragfun")
        {
            if (string(argv[i+1])=="kkp")
                fragfun = new KKP();
            else if (string(argv[i+1])=="pkh")
                fragfun = new PKHFF();
            else if (string(argv[i+1])=="hkns")
                fragfun = new HKNS();
            else if (string(argv[i+1])=="dss")
                fragfun = new DSS(); 
            else
            {
                cerr << "Fragmentation function type " << argv[i+1] << " is not valid!" << endl;
                return -1;
            }
        }
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

    // Default PDF and FF
    if (pdf==NULL)
        pdf = new CTEQ();
    if (fragfun==NULL)
        fragfun = new DSS();
    
    // Read data
    AmplitudeLib N(datafile);

    N.SetFTMethod(ft_method);
    SingleInclusive xs(&N);

    
    
    if (datafile_probe != "")
        datafile_probe = datafile;
    AmplitudeLib N2(datafile_probe);
    if (x0>0)
    {
        N.SetX0(x0); N2.SetX0(x0);
    }
    N2.SetFTMethod(ft_method);

    pdf->SetOrder(order);
    fragfun->SetOrder(order);

    time_t now = time(0);
    string today = ctime(&now);
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# AmplitudeLib v. " << N.Version()  << " running on " << hostname << endl;
    cout <<"# Now is " << today ;
    cout <<"#"<<endl;
	delete[] hostname;

    cout << "# " << N.GetString() << endl;
    cout << "# Fragfun: " << fragfun->GetString() << endl;
    cout << "# PDF: " << pdf->GetString() << endl;
    // Print quarks
    std::vector<Parton> ps;
    ps.push_back(LIGHT);
    //ps.push_back(C);
    ps.push_back(G);
    xs.SetPartons(ps);
    cout <<"# Partons: " ;
    for (unsigned int i=0; i<xs.Partons().size(); i++)
    {
        cout << PartonToString(xs.Partons()[i]) << " ";
    }
    cout << endl;

    //**************** Different operation modes

    if (mode==HYBRID_PT)
    {
        pdf->Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# Single inclusive yield, sqrt(s) = " << sqrts << "GeV" << endl;
        
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << " Producing particle " << ParticleStr(final_particle) <<  endl;
        cout << "# p_T   dN/(d^2 p_T dy) " << endl;
        //cout << "# pt   cteq-partonlevel   ugd-partonlevel " << endl;
        for (double pt=minpt; pt<=maxpt*1.01; pt+=ptstep)
        {
            double result = xs.dHadronMultiplicity_dyd2pt(y, pt, sqrts, fragfun, pdf, final_particle, deuteron);
            cout << pt << " " << result << endl;

        }
    }
    
    else if (mode==HYBRID_PARTON)
    {
		cout <<"# Parton level hybrid formalism" << endl;
		cout <<"# sqrt(s)=" << sqrts << " GeV, y=" << y << endl;
		cout <<"# pdf: " << pdf->GetString() << endl;
		cout <<"# p_T   dN/(d^2 p_T dy)-gluon uquark  dquark  squark  " << endl;
		pdf->Initialize();
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
			double xa = pt*std::exp(-y)/sqrts;
			double ya = std::log(N.X0()/xa);
			double xp =  pt*std::exp(y)/sqrts;
			N.InitializeInterpolation(xa);
				

			double scale = pt;
			double sk = N.S_k(pt, xa, FUNDAMENTAL);
			double sk_adj = N.S_k(pt, xa, ADJOINT);

		
			double partonlevel_u = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, U);
			double partonlevel_d = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, D);
			double partonlevel_s = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, S);
			double partonlevel_g = 1.0/SQR(2.0*M_PI)  * sk_adj * pdf->xq(xp, scale, G);
       
            cout << pt << " " << partonlevel_g << " " << partonlevel_u << " " << partonlevel_d << " " << partonlevel_s  << endl;
        }
		
	}
    
    else if (mode==HYBRID_PT_AVG)
    {
        pdf->Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# <dN / d^2p_T>, sqrt(s) = " << sqrts << "GeV, average over y: " << miny << " - " << maxy  << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   dN" << endl;
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
            double result = 0;
            cerr << "AverageHadronMultiplicity is unimplemented!" << endl;
            //xs.AverageHadronMultiplicity(miny, maxy, pt, sqrts, fragfun, pdf,
            //   deuteron, final_particle);
            cout << pt << " " << result << endl;
        }
    }
    else if (mode==HYBRID_MULTIPLICITY)
    {
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        pdf->Initialize();
        cout << "# Hadron production integrated over pt: " << minpt << " - " << maxpt << endl;
        cout << "# y: " << miny << " - " << maxy << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# sqrt(s)=" << sqrts << " GeV" << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cerr << "HadronMultiplicity is not implemented in v2!" << endl;
        //cout << xs.HadronMultiplicity(miny, maxy, minpt, maxpt, sqrts, fragfun, pdf,
        //    deuteron, final_particle) << endl;
    }
	else if (mode==KTFACT_PT)
    {
		if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# dN/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
		cout << "# Using k_T factorization (gluon production); sigma02 = 1, S_T=1" << endl;
		cout << "# Producing particle " << ParticleStr(final_particle) << endl;
		cout << "# Probe: " << N2.GetString() << endl << "# Target: " << N.GetString() << endl;
		cout << "# NOTICE: results must be multiplied by (\\sigma_0/2)^2/S_T" << endl;
        cout << "# p_T   dN/(d^2 p_T dy) hadronlevel   " << endl;
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
            
            double hadronresult = xs.dHadronMultiplicity_dyd2pt_ktfact(y, pt, sqrts, fragfun, final_particle, &N2);
            cout << pt << " " << hadronresult << endl;// << " " << partonresult << endl;
            
        }
    }
    
    else if (mode==KTFACT_PARTON)
    {
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
		cout << "# Using k_T factorization (gluon production); sigma02 = 1, S_T=1" << endl;
		cout << "# Gluon production (parton level)" << endl;
		cout << "# Probe: " << N2.GetString() << endl << "# Target: " << N.GetString() << endl;
		cout << "# NOTICE: results must be multiplied by (\\sigma_0/2)^2/S_T" << endl;
        cout << "# p_T   dN/(d^2 p_T dy) partonlevel   " << endl;
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
			double partonresult = xs.dHadronMultiplicity_dyd2pt_ktfact_parton(y, pt, sqrts, &N2);
            cout << pt << " " << partonresult << endl;            
        }
    }
    

    return 0;
}
    
    
