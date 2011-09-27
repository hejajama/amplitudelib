/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "../tools/tools.hpp"
#include "../amplitudelib/amplitudelib.hpp"
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
    K,
    SATSCALE,
    GD,
    DSIGMADY,
    PTSPECTRUM,
    F2
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
    string datafile="output.dat";

    if (string(argv[1])=="-help")
    {
        cout << "-y y: set rapidity" << endl;
        cout << "-data datafile" << endl;
        cout << "-x: print amplitude in x space" << endl;
        cout << "-ydep r: print N(r,y) as a function of y" << endl;
        cout << "-k: print amplitude in k space" << endl;
        cout << "-ugd: print unintegrated gluon distribution" << endl;
        cout << "-pt_spectrum: print dN/(d^2 p_T dy)" << endl;
        cout << "-dsigmady: print d\\sigma/dy" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        cout << "-F2 Qsqr" << endl;
        return 0;
    }
    
    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile = argv[i+1];
        else if (string(argv[i])=="-x")
            mode=X;
        else if (string(argv[i])=="-ydep")
        {
            mode=YDEP; r=StrToReal(argv[i+1]);
        }
        else if (string(argv[i])=="-k")
            mode=K;
        else if (string(argv[i])=="-ugd")
            mode=GD;
        else if (string(argv[i])=="-dsigmady")
            mode=DSIGMADY;
        else if (string(argv[i])=="-pt_spectrum")
            mode=PTSPECTRUM;
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
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }

    }

    cout << "# Reading data from file " << datafile << endl;
    AmplitudeLib N(datafile);
    N.InitializeInterpolation(y);
    cout << "# y = " << y << endl;

    if (mode==K)
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
    } else if (mode==X)
    {
        cout <<"# Saturation scale r_s in 1/GeV (N(r_s) = " << Ns <<endl;
        cout <<"### " << N.SaturationScale(y, Ns) << endl;
        cout << "# r [1/GeV]     Amplitude   \\partial_r   \\partial2"
         << " r d ln N / d ln r^2" << endl;
        for (REAL r=N.MinR()*1.01; r<N.MaxR(); r*=1.1)
        {
            cout << r << " " << N.N(r, y) <<  " "
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
        cout <<"# y    r_s [1/GeV]" << endl;
        for (REAL y=0; y < N.MaxY(); y+=0.1)
        {
            cout << y << " " << N.SaturationScale(y, Ns) << endl;
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
        REAL sqrts=7000;
        cout << "#d\\sigma/dy d^2p_T, sqrt(s) = " << sqrts << "GeV" << endl;
        cout << "# p_T   d\\sigma" << endl;
        for (REAL pt=1; pt<6; pt+=0.2)
        {
            REAL result = N.dN_gluon_dyd2pt(pt, y, sqrts);
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
        for (double y=0; y<10; y+=0.1)
        {
            cout << y << " " << Qsqr/(4.0*SQR(M_PI)*ALPHA_e)*
                ( N.ProtonPhotonCrossSection(Qsqr, y, 0)
                + N.ProtonPhotonCrossSection(Qsqr, y, 1) ) << endl;
        }

    }

    else
    {
        cerr << "Unkown mode " << argv[3] << endl;
        return -1;
    }
    return 0;
    

}
