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
#include "../amplitudelib/dis.hpp"
#include "data.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <ctime>
#include <unistd.h>

using namespace Amplitude;
using namespace std;

int main(int argc, char* argv[])
{
    std::stringstream infostr;
    infostr << "# Dipole amplitude analyser " << endl << "# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014 " << endl;
    infostr << "# Command: ";
    for (int i=0; i<argc; i++)
        infostr << argv[i] << " ";
    
    cout << infostr.str() << endl;
    
    gsl_set_error_handler(&ErrHandler);
    
    bool only_charm = false;
	
    
    std::string herafile="";    // if not empty, read Q^2,x,y from datafile

    if ( argc==1 or (string(argv[1])=="-help" or string(argv[1])=="--help")  )
    {
		cout <<"==== USAGE ====" << endl;
        cout << "-xbj: Bjorken x (specify only one)" << endl;
        cout << "-data datafile (bk solution)" << endl;
        cout << "-hera heradatafile [charm,total]" << endl;
        cout << "-x0 x0val: overrides x0 value of the datafile" << endl;
		cout << "-gsl_ft: use GSL to compute Fourier transfer" << endl;
        return 0;
        
    }


    double xbj=-1;
    double x0=-1;
    std::string datafile="";

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-xbj")
            xbj = StrToReal(argv[i+1]);
        else if (string(argv[i])=="-data")
            datafile=argv[i+1];
        else if (string(argv[i])=="-hera")
        {
            herafile=argv[i+1];
            
            if (string(argv[i+2])=="charm")
                only_charm=true;
            else if(string(argv[i+2])=="total")
                only_charm=false;
            else
            {
                cerr << "Unknown option " << argv[i+1] << " for hera datafile" << endl;
                exit(1);
            }
        }
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    if (xbj<0 and herafile=="" )
    {
        cerr << "Set xbj " << LINEINFO << endl;
        return -1;
    }
    
    // Read data
    AmplitudeLib N(datafile, false);
    if (x0>0)
        N.SetX0(x0);

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

    DIS dis(&N);
    
    
    
    if (herafile=="")
    {
        for (double Q2=1; Q2 < 128; Q2*=1.1)
        {
            double light_F2 = dis.F2(Q2, xbj, Amplitude::LIGHT, 0.03);
            double charm_F2 = dis.F2(Q2, xbj, Amplitude::C,1.3528);
            double bottom_F2 = dis.F2(Q2,xbj, Amplitude::B, 4.75);
            cout << xbj << " " << Q2 << " " << light_F2 << " " << charm_F2 << " " << light_F2 + charm_F2 + bottom_F2 << endl;
        }
    }
    else
    {
        Data heradata;
        
        heradata.LoadData(herafile,CHARM); // CHARM flag does not affect anything here
        cout << "Q2 x y s_r" <<endl;
        for (int i=0; i < heradata.NumOfPoints(); i++)
        {
            double Q2 =heradata.Qsqr(i);
            if (Q2 < 1 or Q2 > 1000)
                continue;
            
            
            double x = heradata.xbj(i);
            double y = heradata.y(i);
            double sqrts = std::sqrt( Q2/x);
            
            double sigmar_light=0;
            double sigmar_b=0;
            double sigmar_c =dis.ReducedCrossSection(Q2,x, sqrts, Amplitude::C, 1.3528);
            
            if (only_charm==false)
            {
                sigmar_light =dis.ReducedCrossSection(Q2,x, sqrts, Amplitude::LIGHT, 0.03);
                sigmar_b = dis.ReducedCrossSection(Q2,x, sqrts, Amplitude::B, 4.75);
            }
            cout << Q2 << " " << x << " " << y << " " << sigmar_light+sigmar_c+sigmar_b << endl;
        }
    }
    return 0;
}
    
    
