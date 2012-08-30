/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../tools/tools.hpp"
#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/interpolation.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../fragmentation/kkp.hpp"
#include "../fragmentation/pkhff.hpp"
#include "../fragmentation/hkns.hpp"
#include "../fragmentation/dss.hpp"
#include "../pdf/cteq.hpp"
#include "../amplitudelib/virtual_photon.hpp"
#include "../tools/config.hpp"
#include <iostream>
#include <gsl/gsl_errno.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
using std::string;
const string version = "v. 0.1  2011-xx-xx";

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
    PTSPECTRUM_AVG,
    INT_HADRONPROD,
    F2,
    LOGLOGDER,
    DPS,
    PRINT_FF,
    PRINT_PDF,
    TEST
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "#";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);

    FragmentationFunction *fragfun=NULL;

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
    double minpt=1, maxpt=4;
    Order order=NLO;
    double sqrts=200;
    char dps_mode='a';
    double pt1=0,pt2=0,y1=0,y2=0;
    Hadron final_particle = PI0;    // final state particle in single particle
                                    // production
    Parton parton=U;
    string datafile="amplitude.dat";

    if (string(argv[1])=="-help")
    {
        cout << "-y y: set rapidity" << endl;
        cout << "-xbj bjorken_x (overrides -y)" << endl;
        cout << "-data datafile" << endl;
        cout << "-kspace: data is in k space" << endl;
        cout << "-x: print amplitude (space-indep.)" << endl;
        cout << "-ydep r: print N(r,y) as a function of y" << endl;
        cout << "-x_to_k: FT ampltiude from x to k space" << endl;
        cout << "-k_to_x: FT amplitude from k to x space" << endl;
        cout << "-s_x_to_k: FT S(r)" << endl;
        cout << "-ugd: print unintegrated gluon distribution" << endl;
        cout << "-pt_spectrum p/d pi0/ch/hn: print dN/(d^2 p_T dy), probe is proton or deuteron"
            << " final state pi0/charged/negative hadron" << endl;
        cout << "-pt_spectrum_avg: same as above, but average over y region, must set miny and maxy" << endl;
        cout << "-hadronprod_int p/d pi0/ch/hn: integrated over pt and y range" << endl;
        cout << "-miny y, -maxy y" << endl;
        cout << "-minpt, -maxpt" << endl;
        cout << "-dps p/d pi0/ch/hn a/b/c: doupe parton scattering" << endl;
        cout << "-pt1, -pt2, -y1, -y2: set pt/y for dps calculation" << endl;
        cout << "-dsigmady: print d\\sigma/dy" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        cout << "-F2 Qsqr" << endl;
        cout << "-loglogder: print d ln N / d ln x^2" << endl;
        cout << "-bspline: use bspline interpolation (for noisy data)" << endl;
        cout << "-fragfun [kkp, pkh, hkns, dss]: select fragmentation function" << endl;
        cout << "-sqrts sqrts (in GeV)" << endl;
        cout << "-print_ff [u,d,s,g] [pi0,pim,pip,hm,hp] qsqr" << endl;
        cout << "-print_pdf [u,d,g] qsqr" << endl;
        cout << "-lo: user LO PDF/FF instead of NLO" << endl;
        cout << "-test: run tests" << endl;
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
        else if (string(argv[i])=="-print_pdf")
		{
			mode = PRINT_PDF;
            Qsqr = StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-lo")
        {
			order=LO;
		}
		else if (string(argv[i])=="-test")
			mode=TEST;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }

    cout << "# Reading data from file " << datafile << endl;
    cout <<"# Order: "; if (order==LO) cout << "LO"; else cout << "NLO"; cout << endl;
    CTEQ pdf;
    if (fragfun==NULL)
		fragfun = new DSS(); 
    fragfun->SetOrder(order);
    pdf.SetOrder(order);
    cout << "# PDF: " << pdf.GetString() <<", FF: " << fragfun->GetString() << endl;
    AmplitudeLib N(datafile, kspace);
    N.InitializeInterpolation(y,bspline);
    if (xbj>=0) y = std::log(N.X0()/xbj);
    cout << "# y = " << y << ", x_0 = " << N.X0() << " x = " << N.X0()*std::exp(-y) << endl;

    
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
		double qs = 1.0/N.SaturationScale(y, 0.22);
		cout << "#FT of S(r), Q_s = " << qs << endl;
        double mink = 1e-5; double maxk = 40;// 1.0/N.MinR()*100;
        int kpoints=500;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
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
        cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2 "
         << " r d ln N / d ln r^2" << endl;
        double minr = N.MinR()*1.1; double maxr=N.MaxR()*0.99;
        for (double r=minr; r<maxr; r*=1.03)
        {
            cout << std::scientific << std::setprecision(9) << r << " " << N.N(r, y)  <<  " "
             << N.N(r,y,1) << " " << N.N(r,y,2) <<
             " " << N.LogLogDerivative(r,y) << endl;
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
        cout <<"# Saturation scale N(r_s) = " << Ns << endl;
        cout <<"# y    Q_s [GeV]    d ln Q_s/dy   x" << endl;

        // Solve satscale and save it to array, then interpolate=>get also derivative
        int points = (int)(N.MaxY()/0.1);
        double* rapidities = new double[points];
        double* lnqs = new double[points];

        for (int i=0; i<points; i++)
        {
            rapidities[i]=(double)i*0.1;
            lnqs[i] = std::log(1.0/N.SaturationScale(rapidities[i], Ns));
        }
        Interpolator interp(rapidities, lnqs, points);
        interp.Initialize();

        for (double y=0; y < rapidities[points-1]; y+=0.1)
        {
            cout << y << " " << std::exp(interp.Evaluate(y)) << " "
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
        pdf.Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   dN/(d^2 p_T dy)     parton level yield     F(\\delta)" << endl;
        
        for (double pt=minpt; pt<maxpt; pt+=0.1)
        {
            double result = N.dHadronMultiplicity_dyd2pt(y, pt, sqrts, fragfun, &pdf,
                deuteron, final_particle);
            double partonlevel = N.dHadronMultiplicity_dyd2pt_parton(y, pt, sqrts, &pdf, deuteron);
            double ya = std::log(N.X0() / (pt*std::exp(-y)/sqrts) );
            cout << pt << " " << result << " " << partonlevel << " " << N.S_k(pt, ya)/(4.0*SQR(M_PI)) << endl;
        }
    }
    else if (mode==PTSPECTRUM_AVG)
    {
        pdf.Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# <d\\sigma / d^2p_T>, sqrt(s) = " << sqrts << "GeV, average over y: " << miny << " - " << maxy  << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   d\\sigma" << endl;
        for (double pt=1; pt<4; pt+=0.1)
        {
            double result = N.AverageHadronMultiplicity(miny, maxy, pt, sqrts, fragfun, &pdf,
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
        pdf.Initialize();
        cout << "# Hadron production integrated over pt: " << minpt << " - " << maxpt << endl;
        cout << "# y: " << miny << " - " << maxy << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# sqrt(s)=" << sqrts << " GeV" << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << N.HadronMultiplicity(miny, maxy, minpt, maxpt, sqrts, fragfun, &pdf,
            deuteron, final_particle) << endl;
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
			cout <<"# xf(x_p, u) + xf(x_p, d) = " << pdf.xq(xp, std::max(pt1, pt2), U) + pdf.xq(xp, std::max(pt1, pt2), D) 
			<<", x_p=" << xp << endl;
			
			cout << "# Partonlevel DPS, y1=" << y1 << ", y2=" << y2 << ", pt1=" << pt1 <<", pt2=" << pt2 ;
			double dps_partonlevel = N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, &pdf, deuteron, dps_mode);
			double dps_b =  N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, &pdf, deuteron, 'b');
			double dps_c =  N.DPS_partonlevel(y1, y2, pt1, pt2, sqrts, &pdf, deuteron, 'c');
			cout <<"# partonlevel single inclusive, 1: " << N.dHadronMultiplicity_dyd2pt_parton(y1, pt1, sqrts, &pdf, deuteron) 
				<< " 2: " << N.dHadronMultiplicity_dyd2pt_parton(y2, pt2, sqrts, &pdf, deuteron) << endl;
			cout << "# DPS " << dps_mode <<": " << dps_partonlevel << endl;
			cout <<"# DPS b+c = " << dps_b + dps_c << endl;
        }
        double dps = N.DPSMultiplicity(miny,maxy,1,2,sqrts,fragfun, &pdf, deuteron, final_particle, dps_mode);
        //double single_sqr = N.dHadronMultiplicity_dyd2pt(y1, pt1, sqrts, fragfun, &pdf, deuteron, final_particle)
		//		* N.dHadronMultiplicity_dyd2pt(y2, pt2, sqrts, fragfun, &pdf, deuteron, final_particle);
		cout << "# Hadronlevel DPS " << dps_mode <<": " << endl << dps << endl; //<< " " << single_sqr << " " << dps+single_sqr << endl;
		
		
    }
    
    else if (mode==DSIGMADY)
    {
        double miny=-3;
        double maxy=3;
        int ypoints=30;
        cout << "#d\\sigma/dy, sqrt(s) = 200" << endl;
        cout << "# y     d\\sigma/dy" << endl;
        #pragma omp parallel for
        for (int yind=0; yind<=ypoints; yind++)
        {
            double tmpy = miny + (maxy-miny)/ypoints*yind;
            double result = N.dSigmady_mc(y, 200);
            //double result = N.dSigmadyd2pt(3, 3.0/200.0*std::exp(tmpy), 3.0/200.0*std::exp(-tmpy));
            #pragma omp critical
            {
                cout << tmpy << " " << result << endl;
            }
        }
    }
    else if (mode==F2)
    {
        cout <<"# F_2 at Q^2=" << Qsqr << " GeV^2" << endl;
        VirtualPhoton wf;
        cout << "# Virtual photon wavef params: " << wf.GetParamString() << endl;       
        cout <<"# x   F_2   F_L   scaled_x     y" << endl;
        for(double x=1e-5; x<=N.X0(); x*=1.2)
        {
            // To go smoothly into the photoproduction region, scale
            // x -> x*(1 + 4m_f^2/Qsqr)
            double x2 = x*(1.0+4.0*SQR(0.14)/Qsqr);
            double y = std::log(N.X0()/x2);    // TODO: or x2?
            double xs_l = N.ProtonPhotonCrossSection(Qsqr, y, 0);
            double xs_t = N.ProtonPhotonCrossSection(Qsqr, y, 1);
            cout << x << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t)
                << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l << " "
                << x2 << " " << y << endl;

            
        }
        double xs_l = N.ProtonPhotonCrossSection(Qsqr, 0, 0);
        double xs_t = N.ProtonPhotonCrossSection(Qsqr, 0, 1);
        cout << N.X0() << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t)
                << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l << " "
                << N.X0()*(1.0+4.0*SQR(0.14)/Qsqr) << " " << 0 << endl;

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
		cout <<"# PDF: " << pdf.GetString() << ", Q^2 = " << Qsqr << " GeV^2 " << endl;
		cout << "# x    x*f_u  x*f_d  x*f_s  x_f_g " << endl;
		for (double x=1e-3; x<1; x*=1.05)
		{
			cout << x << " " << pdf.xq(x, q, U) << " " << pdf.xq(x,q,D)
				<< " " << pdf.xq(x,q,S) << " " << pdf.xq(x,q,G) << endl;
		}
	}
	else if (mode==TEST)
	{
		cout << "##### Running tests......" << endl;
		cout <<"### PDF: " << endl;
		pdf.Test();
		cout << endl << "### Fragfun: " << endl;
		fragfun->Test();
	}
    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }

    delete fragfun;
    return 0;
    

}
