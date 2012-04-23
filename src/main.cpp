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
    PRINT_PDF
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
    double sqrts=200;
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
        cout << "-dps: doupe parton scattering, same arguments as -pt_spectrum" << endl;
        cout << "-dsigmady: print d\\sigma/dy" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        cout << "-F2 Qsqr" << endl;
        cout << "-loglogder: print d ln N / d ln x^2" << endl;
        cout << "-bspline: use bspline interpolation (for noisy data)" << endl;
        cout << "-fragfun [kkp, pkh, hkns, dss]: select fragmentation function" << endl;
        cout << "-sqrts sqrts (in GeV)" << endl;
        cout << "-print_ff [u,d,s,g] [pi0,pim,pip,hm,hp] qsqr" << endl;
        cout << "-print_pdf [u,d,g] qsqr" << endl;
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
        }
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
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }

    cout << "# Reading data from file " << datafile << endl;
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
        cout << "# k [GeV]     Amplitude   k/Q_s" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.N_k(tmpk, y);
            #pragma omp critical
            {
                cout <<tmpk << " " << res << " " << tmpk/qs << endl;
            }
        }
    }
    else if (mode==S_X_TO_K)
    {
		double qs = 1.0/N.SaturationScale(y, 0.22);
		cout << "#FT of S(r), Q_s = " << qs << endl;
        double mink = 1e-5; double maxk = 1.0/N.MinR()*100;
        int kpoints=100;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# k [GeV]     Amplitude   k/Q_s" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.S_k(tmpk, y);
            cout <<tmpk << " " << res << " " << tmpk/qs << endl;
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
            #pragma omp critical
            {
                cout <<tmpx << " " << res << endl;
            }
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
        //#pragma omp parallel for schedule(dynamic, 5)
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double result = N.UGD(tmpk, y);
            //#pragma omp critical
            {
                cout << tmpk << " " << result << " " << Alpha_s(SQR(tmpk)) <<endl;
            }
        }
    }

    else if (mode==PTSPECTRUM)
    {
        CTEQ pdf;
        pdf.Initialize();
        if (fragfun==NULL)
        {
            cerr << "Fragfun not spesified!" << endl;
            return -1;
        }
        cout << "# d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   d\\sigma" << endl;
        for (double pt=1; pt<5; pt+=0.1)
        {
            double result = N.dHadronMultiplicity_dyd2pt(y, pt, sqrts, fragfun, &pdf,
                deuteron, final_particle);;
            cout << pt << " " << result << endl;
        }
    }
    else if (mode==PTSPECTRUM_AVG)
    {
        CTEQ pdf;
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
        CTEQ pdf;
        pdf.Initialize();
        double minpt=2; double maxpt=7;
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
        CTEQ pdf; pdf.Initialize();
        double pt1=1, pt2=2, y1=4.2, y2=4.1;
        cout << "# DPS " << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# Fragfun: " << fragfun->GetString() << endl;
        cout << "# sqrt(s)=" << sqrts << " GeV" << endl;
        cout << "# pt1: " << pt1 << " pt2: " << pt2 << " y1: " << y1 << " y2: " << y2 << endl;
        cout <<"# (a)+(c)  (b)   sum" << endl;
        //double dps = N.DPS(y1,y2,pt1,pt2,sqrts, fragfun, &pdf, deuteron, final_particle);
        double dps = N.DPSMultiplicity(2.4,4,1,2,sqrts,fragfun, &pdf, deuteron, final_particle);
        //double single_sqr = N.dHadronMultiplicity_dyd2pt(y1, pt1, sqrts, fragfun, &pdf, deuteron, final_particle)
		//		* N.dHadronMultiplicity_dyd2pt(y2, pt2, sqrts, fragfun, &pdf, deuteron, final_particle);
		cout << dps << endl; //<< " " << single_sqr << " " << dps+single_sqr << endl;
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
		CTEQ pdf;
		pdf.Initialize();
		double q = std::sqrt(Qsqr);
		cout <<"# PDF: " << pdf.GetString() << ", Q^2 = " << Qsqr << " GeV^2 " << endl;
		cout << "# x    x*f_u  x*f_d  x*f_s  x_f_g " << endl;
		for (double x=1e-3; x<1; x*=1.05)
		{
			cout << x << " " << pdf.xq(x, q, U) << " " << pdf.xq(x,q,D)
				<< " " << pdf.xq(x,q,S) << " " << pdf.xq(x,q,G) << endl;
		}
	}
    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }

    delete fragfun;
    return 0;
    

}
