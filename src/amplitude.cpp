/*
 * Dipole amplitude analysis
 * Uses AmplitudeLib to analyse dipole amplitude obtained by solving
 * the BK equation
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2014
 */

#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/config.hpp"
#include "../tools/tools.hpp"
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
    X,
    LOGLOGDER,
    X_TO_K,
    K_TO_X,
    S_X_TO_K,
    SATSCALE
};

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "# Dipole amplitude analyser " << endl << "# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014 " << endl;
    infostr << "# Command: ";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);


    if ( argc==1 or (string(argv[1])=="-help" or string(argv[1])=="--help")  )
    {
		cout <<"==== USAGE ====" << endl;
        cout << "-y or -xbj: set evolution rapidity or Bjorken x (specify only one)" << endl;
        cout << "-data datafile (bk solution)" << endl;
        cout << "-x0 x0val: overrides x0 value of the datafile" << endl;
        cout << "-kspace: data is in k (momentum) space" << endl;
        cout << endl;
        cout << "-x: print amplitude as a function of r or k" << endl;
        cout << "-loglogder: print d ln N / d ln x^2" << endl;
        cout << "-x_to_k: FT N(r)/r^2 from x to k space, compute WW UGD (without prefactors)" << endl;
        cout << "-k_to_x: FT amplitude from k to x space" << endl;
        cout << "-s_x_to_k: FT S(r)" << endl;
        cout << "-satscale Ns, print satscale r_s defined as N(r_s)=Ns" << endl;
        return 0;
        
    }



    double y=-1;
    double xbj=-1;
    double x0=-1;
    double Ns=0.2;
    bool kspace=false;
    std::string datafile="";
    Mode mode=X;

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-y")
            y = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-xbj")
            xbj = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile=argv[i+1];
        else if (string(argv[i])=="-x")
            mode=X;
        else if (string(argv[i])=="-kspace")
            kspace=true;
        else if (string(argv[i])=="-satscale")
        {
            Ns=StrToReal(argv[i+1]);
            mode=SATSCALE;
        }
        else if (string(argv[i])=="-s_x_to_k")
            mode=S_X_TO_K;
        else if (string(argv[i])=="-x_to_k")
            mode=X_TO_K;
        else if (string(argv[i])=="-loglogder")
            mode=LOGLOGDER;
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    if (xbj<0 and y<0 and mode != SATSCALE)
    {
        cerr << "Set xbj or y! " << LINEINFO << endl;
        return -1;
    }
    
    // Read data
    AmplitudeLib N(datafile, kspace);
    if (x0>0)
        N.SetX0(x0);
    if (y>=0)   // User set rapidity
        xbj=N.X0()*std::exp(-y);
    else // User set xbj
        y = std::log(N.X0()/xbj);


    time_t now = time(0);
    string today = ctime(&now);
    
    char *hostname = new char[500];
    gethostname(hostname, 500);
    
    cout <<"#"<<endl<<"# AmplitudeLib v. " << N.Version()  << " running on " << hostname << endl;
    cout <<"# Now is " << today ;
    cout <<"#"<<endl;
	delete[] hostname;

    N.InitializeInterpolation(xbj);
    cout << "# " << N.GetString() << endl;



    //**************** Different operation modes
    if (mode==X)
    {
        cout << "# r [1/GeV]   N(x=" << xbj << " [y=" << y << "])" << endl;
        double minr = N.MinR()*1.1; double maxr=N.MaxR()*0.99;
        for (double r=minr; r<maxr; r*=1.03)
        {
            cout << std::scientific << std::setprecision(9) << r << " " << N.N(r, xbj) << endl;
        }
    }

    else if (mode==S_X_TO_K)
    {
		cout << "#FT of S(r), x=" << xbj << " (y=" << y << ")" << endl;
        double mink = 1e-5; double maxk = 40;
        int kpoints=500;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# k [GeV]     Amplitude    Adj. amplitude" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.S_k(tmpk, xbj, FUNDAMENTAL);
            double adjres = N.S_k(tmpk, xbj, ADJOINT);
            cout <<tmpk << " " << res << " " <<  adjres << endl;
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
            double res = N.N_k_to_x(tmpx, xbj);
            cout <<tmpx << " " << res << endl;
        }
    }

    

    else  if (mode==X_TO_K)
    {
        double sat_n = 0.39346;
		double qs = 1.0/N.SaturationScale(xbj, sat_n);
		cout << "#FT of N(r)/r^2 = WW_UGD" <<  endl;
        double mink = 1e-5; double maxk = 1.0/N.MinR()*100;
        int kpoints=500;
        double kmultiplier = std::pow(maxk/mink, 1.0/(kpoints-1.0));
        cout << "# k [GeV]     WW UGD    FT of S (adjoint)   k/Q_s" << endl;
        for (int kind=0; kind<kpoints; kind++)
        {
            double tmpk = mink*std::pow(kmultiplier, kind);
            double res = N.WW_UGD(tmpk, xbj);
            double ft_s=N.S_k(tmpk, xbj, ADJOINT);
            cout <<tmpk << " " << res << " " << ft_s << " " << tmpk/qs << endl;
        }
    
    }

    else if (mode==SATSCALE)
    {
        cout <<"# Saturation scale N(r^2=2/Qs^2) = " << Ns << endl;
        cout <<"# y    Q_s [GeV]    d ln Q_s^2/dy   x" << endl;

        // Solve satscale and save it to array, then interpolate=>get also derivative
        // Evaluate satscale only at rapidities found from the bk solution file
        // to minimize numerical uncertainty from interpolation

        std::vector<double> rapidities;
        std::vector<double> lnqsqr;

        for (int i=0; i<N.YPoints(); i++)
        {
            double rapidity = N.YValue(i);
            rapidities.push_back(rapidity);
            double x = N.X0()*std::exp(-rapidities[i]);
            double rs = N.SaturationScale(x, Ns);
            lnqsqr.push_back(std::log( 2.0 / SQR(rs) ));    // qs^2
        }
        Interpolator interp(rapidities, lnqsqr);
        interp.Initialize();

        for (int i=0; i < N.YPoints(); i++)
        {
            double y = N.YValue(i);
            cout << y << " " << std::exp(0.5*interp.Evaluate(y)) << " "
                << interp.Derivative(y) <<  " " << N.X0()*std::exp(-y) << endl;
            
        }
    }

    else if (mode==LOGLOGDER)
    {
        cout <<"# d ln N / d ln x^2" << endl;
        cout << "# r [1/GeV]   d ln N / d ln r^2" << endl;
        cout <<"# r    anomalous_dimension" << endl;
        for (double r=N.MinR()*1.1; r<N.MaxR()*0.99; r*=1.03)
        {
            cout << r << " " << N.LogLogDerivative(r, xbj) << endl;
        }

    }

    return 0;
}
    
    
