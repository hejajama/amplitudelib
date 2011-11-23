/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../tools/tools.hpp"
#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/interpolation.hpp"
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
    SATSCALE,
    GD,
    DSIGMADY,
    PTSPECTRUM,
    F2,
    LOGLOGDER
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "#";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);

    Mode mode=X;
    REAL Ns=0.22;
    REAL y=0;
    REAL r=-1;
    REAL Qsqr=10;
    bool kspace=false;
    bool bspline=false;
    bool deuteron=false;
    string datafile="output.dat";

    if (string(argv[1])=="-help")
    {
        cout << "-y y: set rapidity" << endl;
        cout << "-data datafile" << endl;
        cout << "-kspace: data is in k space" << endl;
        cout << "-x: print amplitude (space-indep.)" << endl;
        cout << "-ydep r: print N(r,y) as a function of y" << endl;
        cout << "-x_to_k: FT ampltiudet from x to k space" << endl;
        cout << "-k_to_x: FT amplitude from k to x space" << endl;
        cout << "-ugd: print unintegrated gluon distribution" << endl;
        cout << "-pt_spectrum p/d: print dN/(d^2 p_T dy), probe is proton or deuteron" << endl;
        cout << "-dsigmady: print d\\sigma/dy" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        cout << "-F2 Qsqr" << endl;
        cout << "-loglogder: print d ln N / d ln x^2" << endl;
        cout << "-bspline: use bspline interpolation (for noisy data)" << endl;
        return 0;
    }
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile = argv[i+1];
        else if (string(argv[i])=="-kspace")
            kspace=true;
        else if (string(argv[i])=="-x")
            mode=X;
        else if (string(argv[i])=="-ydep")
        {
            mode=YDEP; r=StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-x_to_k")
            mode=X_TO_K;
        else if (string(argv[i])=="-k_to_x")
            mode=K_TO_X;
        else if (string(argv[i])=="-ugd")
            mode=GD;
        else if (string(argv[i])=="-dsigmady")
            mode=DSIGMADY;
        else if (string(argv[i])=="-pt_spectrum")
        {
            mode=PTSPECTRUM;
            if (string(argv[i+1])=="p")
                deuteron=false;
            else if (string(argv[i+1])=="d")
                deuteron=true;
            else
            {
                cerr << "Invalid probe particle type " << argv[i+1] << endl;
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
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }

    cout << "# Reading data from file " << datafile << endl;
    AmplitudeLib N(datafile, kspace);
    N.InitializeInterpolation(y,bspline);
    cout << "# y = " << y << endl;

    if (mode==X_TO_K)
    {
        REAL mink = 1e-5; REAL maxk = 1.0/N.MinR()*100;
        int kpoints=100;
        REAL kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# k [GeV]     Amplitude  " << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            REAL tmpk = mink*std::pow(kmultiplier, kind);
            REAL res = N.N_k(tmpk, y);
            #pragma omp critical
            {
                cout <<tmpk << " " << res << endl;
            }
        }
    }
    else if (mode==K_TO_X)
    {
        REAL minx = 1e-6; REAL maxx = 50;
        int xpoints=100;
        REAL xmultiplier = std::pow(maxx/minx, 1.0/(xpoints-1.0));
        cout << "# x [1/GeV]     Amplitude  " << endl;
        for (int xind=0; xind<xpoints; xind++)
        {
            REAL tmpx = minx*std::pow(xmultiplier, xind);
            REAL res = N.N_k_to_x(tmpx, y);
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
        cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2"
         << " r d ln N / d ln r^2" << endl;
        REAL minr = N.MinR()*1.1; REAL maxr=N.MaxR()*0.99;
        for (REAL r=minr; r<maxr; r*=1.03)
        {
            cout << std::scientific << std::setprecision(7) << r << " " << N.N(r, y)  <<  " "
             << N.N(r,y,1) << " " << N.N(r,y,2) <<
             " " << N.LogLogDerivative(r,y) << endl;
        }
    }
    else if (mode==YDEP)
    {
        cout << "# N(r=" << r <<", y) as sa function of y" << endl;
        for (REAL y=0; y<N.MaxY(); y+=0.1)
        {
            cout << y << " " << N.N(r,y) << endl;
        }

    }
    else if (mode==SATSCALE)
    {
        cout <<"# Saturation scale N(r_s) = " << Ns << endl;
        cout <<"# y    Q_s [GeV]    d ln Q_s/dy" << endl;

        // Solve satscale and save it to array, then interpolate=>get also derivative
        int points = (int)(N.MaxY()/0.1);
        REAL* rapidities = new REAL[points];
        REAL* lnqs = new REAL[points];

        for (int i=0; i<points; i++)
        {
            rapidities[i]=(REAL)i*0.1;
            lnqs[i] = std::log(1.0/N.SaturationScale(rapidities[i], Ns));
        }
        Interpolator interp(rapidities, lnqs, points);
        interp.Initialize();

        for (REAL y=0; y < rapidities[points-1]; y+=0.1)
        {
            cout << y << " " << std::exp(interp.Evaluate(y)) << " "
                << interp.Derivative(y) << endl;
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
        REAL mink=0.1; REAL maxk=20;
        int kpoints=250;
        REAL kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# UGD" << endl << "# k_T [GeV]   UGD   \\alpha_s(k)" << endl;
        //#pragma omp parallel for schedule(dynamic, 5)
        for (int kind=0; kind<kpoints; kind++)
        {
            REAL tmpk = mink*std::pow(kmultiplier, kind);
            REAL result = N.UGD(tmpk, y);
            //#pragma omp critical
            {
                cout << tmpk << " " << result << " " << Alpha_s(SQR(tmpk)) <<endl;
            }
        }
    }

    else if (mode==PTSPECTRUM)
    {
        REAL sqrts=200;
        cout << "#d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
        cout << "# Probe: "; if (deuteron) cout <<"deuteron"; else cout <<"proton"; cout << endl;
        cout << "# p_T   d\\sigma" << endl;
        for (REAL pt=0.5; pt<5; pt+=0.1)
        {
            REAL result = N.dHadronMultiplicity_dyd2pt(y, pt, sqrts, deuteron);
            cout << pt << " " << result << endl;
        }
    }
    
    else if (mode==DSIGMADY)
    {
        REAL miny=-3;
        REAL maxy=3;
        int ypoints=30;
        cout << "#d\\sigma/dy, sqrt(s) = 200" << endl;
        cout << "# y     d\\sigma/dy" << endl;
        #pragma omp parallel for
        for (int yind=0; yind<=ypoints; yind++)
        {
            REAL tmpy = miny + (maxy-miny)/ypoints*yind;
            REAL result = N.dSigmady_mc(y, 200);
            //REAL result = N.dSigmadyd2pt(3, 3.0/200.0*std::exp(tmpy), 3.0/200.0*std::exp(-tmpy));
            #pragma omp critical
            {
                cout << tmpy << " " << result << endl;
            }
        }
    }
    else if (mode==F2)
    {
        cout <<"# F_2 at Q^2=" << Qsqr << " GeV^2" << endl;
        cout <<"# x   F_2   F_L   scaled_x     y" << endl;
        for (double x=1e-6; x<=1e-2; x*=1.1)
        {
            // To go smoothly into the photoproduction region, scale
            // x -> x*(1 + 4m_f^2/Qsqr)
            double x2 = x*(1.0+4.0*SQR(0.14)/Qsqr);
            double y = std::log(N.X0()/x);    // TODO: or x2?
            double xs_l = N.ProtonPhotonCrossSection(Qsqr, y, 0);
            double xs_t = N.ProtonPhotonCrossSection(Qsqr, y, 1);
            cout << x << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*(xs_l+xs_t)
                << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*xs_l << " "
                << x2 << " " << y << endl;
        }

    }

    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }
    return 0;
    

}
