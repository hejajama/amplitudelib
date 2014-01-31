/*
 * AmplitudeLib, reads output of the BK equation solver 
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
 */

#include "../tools/tools.hpp"
#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/interpolation.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../fragmentation/kkp.hpp"
#include "../fragmentation/pkhff.hpp"
#include "../fragmentation/hkns.hpp"
#include "../fragmentation/dss.hpp"
#include "../pdf/ugdpdf.hpp"
#include "../pdf/cteq.hpp"
#include "../pdf/mrst.hpp"
#include "../pdf/eps09.hpp"
#include "../amplitudelib/virtual_photon.hpp"
#include "../tools/config.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <unistd.h>

using namespace std;


enum Mode
{
    X,
    YDEP,
    X_TO_K,
    K_TO_X,
    S_X_TO_K,
    SATSCALE,
    GD,
    DSIGMADY,
    PTSPECTRUM,
    PTSPECTRUM_PARTON,
    PTSPECTRUM_AVG,
    PTSPECTRUM_KTFACT,
    PTSPECTRUM_KTFACT_PARTON,
    INT_HADRONPROD,
    F2,
    LOGLOGDER,
    DPS,
    PRINT_FF,
    PRINT_PDF,
    PRINT_PDF_Q,
    TEST
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "# amplitudeLib   (c) Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013 " << endl;
    infostr << "# Command: ";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);

    FragmentationFunction *fragfun=NULL;
    PDF *pdf=NULL;

    Mode mode=X;
    double Ns=0.22;
    double y=0;
    double xbj=-1;
    double r=-1;
    double Qsqr=10;
    bool kspace=false;
    bool bspline=false;
    bool deuteron=false;
    double miny=3; double maxy=4;
    double minpt=1, maxpt=8;
    double ptstep = 0.1;
    Order order=NLO;
    double sqrts=200;
    char dps_mode='a';
    double pt1=0,pt2=0,y1=0,y2=0;
    double sigma02=1.0;
    RUNNING_ALPHAS as = RUNNING;
    FT_METHOD ft = ACC_SERIES;
    Hadron final_particle = PI0;    // final state particle in single particle
                                    // production
    Parton parton=U;
    string datafile="amplitude.dat";
    string ktfact_datafile2="";
    double x0=-1;	// use default

    if (string(argv[1])=="-help")
    {
		cout <<"==== General parameters ====" << endl;
        cout << "-y y: set rapidity" << endl;
        cout << "-xbj bjorken_x (overrides -y)" << endl;
        cout << "-data datafile (bk solution)" << endl;
        cout << "-x0 x0val: overrides x0 value of the datafile" << endl;
        cout << "-gsl_ft: use GSL numerical integral to compute FT for S_k " << endl;
        cout << "-bessel_ft: use in principle more advanced accelerated summation method to compute FT for S_k (default)" << endl;
        cout << "-kspace: data is in k (momentum) space" << endl;
        cout << "-bspline: use bspline interpolation (for noisy data) [EXPERIMENTAL]" << endl;
        cout << "-loglogder: print d ln N / d ln x^2" << endl;
        cout << "-sqrts sqrts (in GeV)" << endl;
        
        cout << endl << "==== Study dipole amplitude " << endl;
        cout << "-x: print amplitude as a function of r or k" << endl;
        cout << "-ydep r: print N(r,y) as a function of y" << endl;
        cout << "-x_to_k: FT N(r)/r^2 from x to k space" << endl;
        cout << "-k_to_x: FT amplitude from k to x space" << endl;
        cout << "-s_x_to_k: FT S(r)" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        
        cout << endl <<"==== Single inclusive ====" << endl;
        cout << "-pt_spectrum p/d pi0/ch/hn: print dN/(d^2 p_T dy), probe is proton or deuteron"
            << " final state pi0/charged/negative hadron" << endl;
        cout << "-pt_spectrum_parton: dN/(d^2 p_t dY) for quark/gluon scattering" << endl;
        cout << "-pt_spectrum_avg: same as above, but average over y region, must set miny and maxy" << endl;
        cout << "-pt_spectrum_ktfact [final particle]: hadron production dN/(d^2 p_T dy) using k_T factorization" << endl;
        cout << "-pt_spectrum_ktfact_parton: dN/(d^2 p_t dY) for gluon production" << endl;
        cout << "-ktfact_probe datafile: in asymmetric collisions dipole amplitude for probe " << endl;
        cout << "-fixed_alphas: use fixed coupling in ktfactorization" << endl;
        cout << "-sigma02 val: proton or target area; ktfactorization results are multiplied by this" << endl;
        cout << "-dsigmady: print d\\sigma/dy" << endl;
        cout << "-hadronprod_int p/d pi0/ch/hn: integrated over pt and y range" << endl;
        cout << "-miny y, -maxy y" << endl;
        cout << "-minpt, -maxpt, -ptstep" << endl;
        
        cout << endl << "==== Double inclusive ====" << endl;
        cout << "-dps p/d pi0/ch/hn a/b/c: double parton scattering" << endl;
        cout << "-pt1, -pt2, -y1, -y2: set pt/y for dps calculation" << endl;
        
        
        cout << endl << "==== PDF and FF ====" << endl;
        cout << "-pdf [ctreq, ugd, ugd_fixed, eps09] [params]   , parameters for ugdpdf: amplitudefile sigma0/2 (ugd_fixed is fixed alphas); for eps09: A" << endl;
        cout << "-lo: use LO PDF and FF;  -nlo: use NLO PDF and FF (default) "<< endl;
        cout << "-fragfun [kkp, pkh, hkns, dss]: select fragmentation function" << endl;
        cout << "-ugd: print unintegrated gluon distribution computed from dipole amplitude" << endl;
        cout << "-print_ff [u,d,s,g] [pi0,pim,pip,hm,hp] qsqr" << endl;
        cout << "-print_pdf qsqr" << endl;
        cout << "-print_pdf_q x" << endl;
        
        
        cout << endl << "==== Misc ==== " << endl;
        cout << "-F2 Qsqr   compute structure function F2 and reduced cross section" << endl;
        cout << "-test: run tests" << endl;
        cout << endl << "==== Notes: ==== " << endl;
        cout << "* All dimensionfull values are GeV^n" << endl;
        return 0;
    }

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xbj")
            xbj = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-miny")
            miny = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxy")
            maxy = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-minpt")
            minpt = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ptstep")
			ptstep = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-maxpt")
            maxpt = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile = argv[i+1];
        else if (string(argv[i])=="-kspace")
            kspace=true;
        else if (string(argv[i])=="-x")
            mode=X;
        else if (string(argv[i])=="-sqrts")
            sqrts = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-ydep")
        {
            mode=YDEP; r=StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-x_to_k")
            mode=X_TO_K;
        else if (string(argv[i])=="-k_to_x")
            mode=K_TO_X;
        else if (string(argv[i])=="-s_x_to_k")
			mode=S_X_TO_K;
		else if (string(argv[i])=="-gsl_ft")
			ft = GSL;
		else if (string(argv[i])=="-bessel_ft")
			ft = ACC_SERIES;
        else if (string(argv[i])=="-ugd")
            mode=GD;
        else if (string(argv[i])=="-dsigmady")
            mode=DSIGMADY;
        else if (string(argv[i])=="-pt_spectrum" or string(argv[i])=="-pt_spectrum_avg")
        {
            if (string(argv[i])=="-pt_spectrum")
                mode=PTSPECTRUM;
            else
                mode=PTSPECTRUM_AVG;
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
			mode=PTSPECTRUM_PARTON;
        else if (string(argv[i])=="-pt_spectrum_ktfact")
		{
			mode = PTSPECTRUM_KTFACT;
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
			mode=PTSPECTRUM_KTFACT_PARTON;
		}
        
        else if (string(argv[i])=="-ktfact_probe")
        {
			ktfact_datafile2=argv[i+1];
		}
			
        else if (string(argv[i])=="-hadronprod_int")
        {
            mode=INT_HADRONPROD;
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
        else if (string(argv[i])=="-dps")
        {
            mode=DPS;
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
            dps_mode = string(argv[i+3])[0];
        }
        else if (string(argv[i])=="-pt1")
			pt1 = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-pt2")
			pt2 = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-y1")
			y1 = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-y2")
			y2 = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-bspline")
            bspline=true;
        else if (string(argv[i])=="-satscale")
        {
            mode=SATSCALE;
            Ns = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-F2")
        {
            mode=F2;
            Qsqr = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-loglogder")
            mode = LOGLOGDER;
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
			else if (string(argv[i+1])=="ugd" or string(argv[i+1])=="ugd_fixed")
			{
				AmplitudeLib* ugdn = new AmplitudeLib(string(argv[i+2]));
				pdf = new UGDPDF(ugdn, StrToReal(argv[i+3]));
				if (string(argv[i+1])=="ugd_fixed")
					ugdn->SetRunningCoupling(FIXED);
			}
			else if (string(argv[i+1])=="eps09")
			{
				pdf = new EPS09();
				pdf->SetA(StrToInt(argv[i+2]));
			}
			else
			{
				cerr << "Unknown PDF type " << argv[i+1] << endl;
				return -1;
			}
		}
        else if (string(argv[i])=="-print_ff")
        {
            mode =  PRINT_FF;
            if (string(argv[i+1])=="u")
                parton=U;
            else if (string(argv[i+1])=="d")
                parton=D;
            else if (string(argv[i+1])=="s")
                parton=S;
            else if (string(argv[i+1])=="g")
                parton=G;
            else
            {
                cerr << "Parton " << string(argv[i+1]) << " is unknown! " << endl;
                return -1;
            }
            if (string(argv[i+2])=="pi0")
                final_particle = PI0;
            else if (string(argv[i+2])=="pim")
                final_particle = PIM;
            else if (string(argv[i+2])=="pip")
                final_particle = PIP;
            else if (string(argv[i+2])=="hm")
                final_particle = HM;
            else if (string(argv[i+2])=="hp")
                final_particle = HP;
            else
            {
                cerr << "Hadron " << string(argv[i+1]) << " is not supported" << endl;
                return -1;
            }

            Qsqr = StrToReal(string(argv[i+3]));
        }
        else if (string(argv[i])=="-print_pdf" )
		{
			mode = PRINT_PDF;
            Qsqr = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-print_pdf_q")
		{
			mode = PRINT_PDF_Q;
			xbj = StrToReal(argv[i+1]);
		}
        else if (string(argv[i])=="-lo")
			order=LO;
		else if (string(argv[i])=="-nlo")
			order=NLO;
		else if (string(argv[i])=="-sigma02")
			sigma02 = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-test")
			mode=TEST;
		else if (string(argv[i])=="-x0")
			x0 = StrToReal(argv[i+1]);
		else if (string(argv[i])=="-fixed_alphas")
			as = FIXED;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }

    cout << "# Reading data from file " << datafile << endl;
   
    AmplitudeLib N(datafile, kspace);
    N.SetFTMethod(ft);
    
    time_t now = time(0);
    string today = ctime(&now);
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# AmplitudeLib v. " << N.Version()  << " running on " << hostname << endl;
    cout <<"# Now is " << today ;
    cout <<"#"<<endl;
	delete[] hostname;
	
	cout <<"# FT method: "; if (N.GetFTMethod()==GSL) cout << "GSL"; else if (N.GetFTMethod()==ACC_SERIES) cout << "Acc Series"; else cout << "Unknown"; cout << endl;
    
   
    cout <<"# Order: "; if (order==LO) cout << "LO"; else cout << "NLO"; cout << endl;
	if (pdf==NULL)
		pdf = new CTEQ();
    if (fragfun==NULL)
		fragfun = new DSS(); 
    fragfun->SetOrder(order);
    pdf->SetOrder(order);
    
    cout << "# PDF: " << pdf->GetString() <<", FF: " << fragfun->GetString() << endl;
   
    if ((mode==PTSPECTRUM_KTFACT or mode==PTSPECTRUM_KTFACT_PARTON) and ktfact_datafile2 != "")
		datafile=ktfact_datafile2;
	
    AmplitudeLib N2(datafile);
    N.SetSigma02(sigma02); N2.SetSigma02(sigma02);
    N2.SetFTMethod(ft);
    N.SetRunningCoupling(as); N2.SetRunningCoupling(as);
    if (x0>0) { N.SetX0(x0); N2.SetX0(x0); }
    N.InitializeInterpolation(y,bspline);
    if (xbj>=0) y = std::log(N.X0()/xbj);
    if (N.GetRunningCoupling()==FIXED)
		cout << "# Fixed alphas = " << N.Alphas(1) << endl;
	else
		cout << "# Running alphas" << endl;
    cout << "# y = " << y << ", x_0 = " << N.X0() << endl;

	
    
    if (mode==X_TO_K)
    {
		double qs = 1.0/N.SaturationScale(y, 0.22);
		cout << "#FT of N(r)/r^2, Q_s = " << qs << endl;
        double mink = 1e-5; double maxk = 1.0/N.MinR()*100;
        int kpoints=100;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# k [GeV]     Amplitude N(k)    FT of S   k/Q_s" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.N_k(tmpk, y);
            double ft_s=N.S_k(tmpk, y);
            cout <<tmpk << " " << res << " " << ft_s << " " << tmpk/qs << endl;
        }
    }
    else if (mode==S_X_TO_K)
    {
		double qs = std::sqrt(2.0)/N.SaturationScale(y, 0.393469);
		cout << "#FT of S(r), Q_s = " << qs << endl;
        double mink = 1e-5; double maxk = 40;// 1.0/N.MinR()*100;
        int kpoints=500;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        N.InitializeInterpolation(y);
        cout << "# k [GeV]     Amplitude    Adj. amplitude    k/Q_s  " << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.S_k(tmpk, y);
            double adjres = N.S_k(tmpk, y, true);
            cout <<tmpk << " " << res << " " <<  adjres << " " << tmpk/qs << endl;
        }
    }
    else if (mode==K_TO_X)
    {
        double minx = 1e-6; double maxx = 50;
        int xpoints=100;
        double xmultiplier = std::pow(maxx/minx, 1.0/(xpoints-1.0));
        cout << "# x [1/GeV]     Amplitude  " << endl;
        for (int xind=0; xind<xpoints; xind++)
        {
            double tmpx = minx*std::pow(xmultiplier, xind);
            double res = N.N_k_to_x(tmpx, y);
            cout <<tmpx << " " << res << endl;
        }
    }

    else if (mode==X)
    {
        cout <<"# Saturation scale r_s in 1/GeV / k_s in GeV (N(r_s) = " << Ns <<endl;
        cout <<"### " << N.SaturationScale(y, Ns) << endl;
        //cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2 "
        // << " r d ln N / d ln r^2" << endl;
        cout << "# r   N" << endl;
        double minr = N.MinR()*1.1; double maxr=N.MaxR()*0.99;
        for (double r=minr; r<maxr; r*=1.03)
        {
            cout << std::scientific << std::setprecision(9) << r << " " << N.N(r, y) << endl;/* << " "
             << N.N(r,y,1) << " " << N.N(r,y,2) <<
             " " << N.LogLogDerivative(r,y) << endl;*/
        }
    }
    else if (mode==YDEP)
    {
        cout << "# N(r=" << r <<", y) as sa function of y" << endl;
        for (double y=0; y<N.MaxY(); y+=0.1)
        {
            cout << y << " " << N.N(r,y) << endl;
        }

    }
    else if (mode==SATSCALE)
    {
        cout <<"# Saturation scale N(r^2=2/Qs^2) = " << Ns << endl;
        cout <<"# y    Q_s [GeV]    d ln Q_s/dy   x" << endl;

        // Solve satscale and save it to array, then interpolate=>get also derivative
        double ystep=0.1;
        int points = (int)(N.MaxY()/ystep);
        double* rapidities = new double[points];
        double* lnqs = new double[points];

        for (int i=0; i<points; i++)
        {
            rapidities[i]=(double)i*ystep;
            double rs = N.SaturationScale(rapidities[i], Ns);
            lnqs[i] = std::log( 2.0 / SQR(rs) );    // qs^2
            //lnqs[i] = std::log(1.0/N.SaturationScale(rapidities[i], Ns));
        }
        Interpolator interp(rapidities, lnqs, points);
        interp.Initialize();

        for (double y=0; y < rapidities[points-1]; y+=ystep)
        {
            cout << y << " " << std::exp(0.5*interp.Evaluate(y)) << " "
                << interp.Derivative(y) <<  " " << N.X0()*std::exp(-y) << endl;
            
        }

        delete[] rapidities;
        delete[] lnqs;
    }

    else if (mode==LOGLOGDER)
    {
        cout <<"# d ln N / d ln x^2" << endl;
        cout <<"# Saturation scale r_s in 1/GeV / k_s in GeV (N(r_s) = " << Ns <<endl;
        cout <<"### " << N.SaturationScale(y, Ns) << endl;
        cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2"
         << " r d ln N / d ln r^2" << endl;
        cout <<"# r    anomalous_dimension" << endl;
        for (double x=N.MinR()*1.1; x<N.MaxR()*0.99; x*=1.03)
        {
            cout << x << " " << N.LogLogDerivative(x, y) << endl;
        }

    }
    
    else if (mode==GD)
    {
        double mink=0.1; double maxk=20;
        int kpoints=250;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# UGD" << endl << "# k_T [GeV]   UGD   \\alpha_s(k)" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double result = N.UGD(tmpk, y);
                cout << tmpk << " " << result << " " << Alpha_s(SQR(tmpk)) <<endl;
        }
    }

    else if (mode==PTSPECTRUM)
    {
        pdf->Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << " Producing particle " << ParticleStr(final_particle) <<  endl;
        cout << "# p_T   dN/(d^2 p_T dy) " << endl;
        //cout << "# pt   cteq-partonlevel   ugd-partonlevel " << endl;
        for (double pt=minpt; pt<=maxpt*1.01; pt+=ptstep)
        {
            double result = N.dHadronMultiplicity_dyd2pt(y, pt, sqrts, fragfun, pdf,deuteron, final_particle);
            cout << pt << " " << result << endl;

        }
    }
    
    else if (mode==PTSPECTRUM_PARTON)
    {
		cout <<"# Parton level hybrid formalism" << endl;
		cout <<"# sqrt(s)=" << sqrts << " GeV" << endl;
		cout <<"# pdf: " << pdf->GetString() << endl;
		cout <<"# p_T   dN/(d^2 p_T dy)-gluon uquark  dquark  squark  " << endl;
		pdf->Initialize();
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
			double xa = pt*std::exp(-y)/sqrts;
			double ya = std::log(N.X0()/xa);
			double xp =  pt*std::exp(y)/sqrts;
			N.InitializeInterpolation(ya);
				

			double scale = pt;
			double sk = N.S_k(pt, ya);
			double sk_adj = N.S_k(pt, ya, true);

		
			double partonlevel_u = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, U);
			double partonlevel_d = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, D);
			double partonlevel_s = 1.0/SQR(2.0*M_PI)  * sk * pdf->xq(xp, scale, S);
			double partonlevel_g = 1.0/SQR(2.0*M_PI)  * sk_adj * pdf->xq(xp, scale, G);
       
            cout << pt << " " << partonlevel_g << " " << partonlevel_u << " " << partonlevel_d << " " << partonlevel_s  << endl;
        }
		
	}
    
    else if (mode==PTSPECTRUM_AVG)
    {
        pdf->Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# <d\\sigma / d^2p_T>, sqrt(s) = " << sqrts << "GeV, average over y: " << miny << " - " << maxy  << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   d\\sigma" << endl;
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
            double result = N.AverageHadronMultiplicity(miny, maxy, pt, sqrts, fragfun, pdf,
                deuteron, final_particle);
            cout << pt << " " << result << endl;
        }
    }
    else if (mode==INT_HADRONPROD)
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
        cout << N.HadronMultiplicity(miny, maxy, minpt, maxpt, sqrts, fragfun, pdf,
            deuteron, final_particle) << endl;
    }
	else if (mode==PTSPECTRUM_KTFACT)
    {
		if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
		cout << "# Using k_T factorization (gluon production); sigma02 = " << N.Sigma02() << " GeV^-2" << endl;
		cout << "# Producing particle " << ParticleStr(final_particle) << endl;
		cout << "# Probe: " << N2.GetString() << endl << "# Target: " << N.GetString() << endl;
		cout << "# NOTICE: results must be multiplied by (\\sigma_0/2)^2/S_T" << endl;
        cout << "# p_T   dN/(d^2 p_T dy) hadronlevel   " << endl;
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
            
            double hadronresult = N.dHadronMultiplicity_dyd2pt_ktfact(y, pt, sqrts, fragfun, final_particle, &N2);
            cout << pt << " " << hadronresult << endl;// << " " << partonresult << endl;
            
        }
    }
    
    else if (mode==PTSPECTRUM_KTFACT_PARTON)
    {
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
		cout << "# Using k_T factorization (gluon production); sigma02 = " << N.Sigma02() << " GeV^-2" << endl;
		cout << "# Gluon production (parton level)" << endl;
		cout << "# Probe: " << N2.GetString() << endl << "# Target: " << N.GetString() << endl;
		cout << "# NOTICE: results must be multiplied by (\\sigma_0/2)^2/S_T" << endl;
        cout << "# p_T   dN/(d^2 p_T dy) partonlevel   " << endl;
        
        for (double pt=minpt; pt<=maxpt; pt+=ptstep)
        {
			double partonresult = N.dHadronMultiplicity_dyd2pt_ktfact_parton(y, pt, sqrts, &N2);
            cout << pt << " " << partonresult << endl;            
        }
    }
    
    else if (mode==DPS)
    {
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# DPS " << dps_mode << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# sqrt(s)=" << sqrts << " GeV" << endl;
        cout << "# Deuteron: "<< deuteron << endl;

        
        if (pt1>0 and pt2>0 and y1>0 and y2>0)
        {
			double xp = (pt1*std::exp(y1) + pt2*std::exp(y2))/sqrts;
			double xa = (pt1*std::exp(-y1) + pt2*std::exp(-y2))/sqrts;
			cout << "# x_p=" << xp <<", x_a=" << xa << endl;
			cout <<"# Q_s = " << 1.0/N.SaturationScale(std::log(N.X0()/xa), 0.22) << " GeV" << endl;
			cout <<"# xf(x_p, u) + xf(x_p, d) = " << pdf->xq(xp, std::max(pt1, pt2), U) + pdf->xq(xp, std::max(pt1, pt2), D) 
			<<", x_p=" << xp << endl;
			
			cout << "# Partonlevel DPS, y1=" << y1 << ", y2=" << y2 << ", pt1=" << pt1 <<", pt2=" << pt2 << endl;
			double dps_partonlevel = N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, pdf, deuteron, dps_mode);
			double dps_b =  N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, pdf, deuteron, 'b');
			double dps_c =  N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, pdf, deuteron, 'c');
			cout <<"# partonlevel single inclusive, 1: " << N.dHadronMultiplicity_dyd2pt_parton(y1, pt1, sqrts, pdf, deuteron) 
				<< " 2: " << N.dHadronMultiplicity_dyd2pt_parton(y2, pt2, sqrts, pdf, deuteron) << endl;
			cout << "# DPS " << dps_mode <<": " << dps_partonlevel << endl;
			cout <<"# DPS b+c = " << dps_b + dps_c << endl;
        }
        double dps = N.DPSMultiplicity(miny,maxy,1,2,sqrts,fragfun, pdf, deuteron, final_particle, dps_mode);
        //double single_sqr = N.dHadronMultiplicity_dyd2pt(y1, pt1, sqrts, fragfun, pdf, deuteron, final_particle)
		//		* N.dHadronMultiplicity_dyd2pt(y2, pt2, sqrts, fragfun, pdf, deuteron, final_particle);
		cout << "# Hadronlevel DPS " << dps_mode <<": " << endl << dps << endl; //<< " " << single_sqr << " " << dps+single_sqr << endl;
		
		
    }
    
    else if (mode==DSIGMADY)
    {
        cerr << "Not implemented! "  << LINEINFO << endl;
    }
    else if (mode==F2)
    {
        cout <<"# F_2 at Q^2=" << Qsqr << " GeV^2" << " sqrts " << sqrts << " GeV" << endl;
        VirtualPhoton wf;
        cout << "# Virtual photon wavef params: " << wf.GetParamString() << endl;       
        cout <<"# x   F_2   F_L  reduced_xs  scaled_x     y" << endl;
        for(double x=1e-5; x<=N.X0(); x*=1.2)
        {
            // To go smoothly into the photoproduction region, scale
            // x -> x*(1 + 4m_f^2/Qsqr)
            double x2 = x*(1.0+4.0*SQR(0.14)/Qsqr);
            double evol_y = std::log(N.X0()/x2);
            
            cout << x << " " << N.F2(Qsqr, evol_y) << " " << N.FL(Qsqr, evol_y)
				<< " " << N.ReducedCrossSection(Qsqr, evol_y, sqrts) << " " << x2 << " " << evol_y << endl;       
        }

    }

    else if (mode == PRINT_FF)
    {
        cout << "# Fragmentation function " << fragfun->GetString() << " at "
         << "Q^2 = " << Qsqr << " GeV^2" << endl;
         cout << "# x  D   x*D" << endl;

        for (double x=0.1; x<1; x*=1.05)
        {
            cout << x << " " << fragfun->Evaluate(parton, final_particle, x, std::sqrt(Qsqr))
             << " " << x*fragfun->Evaluate(parton, final_particle, x, std::sqrt(Qsqr)) << endl;
            
        }         

    }
    
    else if (mode == PRINT_PDF)
    {
		double q = std::sqrt(Qsqr);
		cout <<"# PDF: " << pdf->GetString() << ", Q^2 = " << Qsqr << " GeV^2 " << endl;
		cout << "# x    x*f_u  x*f_d  x*f_s  x*f_g " << endl;
		for (double x=1e-5; x<1; x*=1.05)
		{
			cout << x << " " << pdf->xq(x, q, U) << " " << pdf->xq(x,q,D)
				<< " " << pdf->xq(x,q,S) << " " << pdf->xq(x,q,G) << endl;
		}
	}
	
	else if (mode == PRINT_PDF_Q)
	{
		cout <<"# PDF: " << pdf->GetString() <<", x=" << xbj << endl;
		cout << "# Q^2 [GeV^2]   x*f_u   x*f_d   x*f_s    x*f_g" << endl;
		for (double q=pdf->MinQ(); q*q<std::min(1000.0,pdf->MaxQ()); q*=1.05)
		{
			cout << q*q << " " << pdf->xq(xbj, q, U) << " " << pdf->xq(xbj,q,D)
				<< " " << pdf->xq(xbj,q,S) << " " << pdf->xq(xbj,q,G) << endl;
		}
	}
	

	else if (mode==TEST)
	{
		cout << "##### Running tests......" << endl;
		cout <<"### PDF: " << endl;
		pdf->Test();
		cout << endl << "### Fragfun: " << endl;
		fragfun->Test();
	}
    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }

	if (fragfun != NULL)
		delete fragfun;
    if (pdf != NULL)
		delete pdf;
    return 0;
    

}
