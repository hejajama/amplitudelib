/*
 * Single inclusive cross section calculations using the
 * AmplitudeLib
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
        cout << "-y rapidity";
        cout << "-sqrts center-of-mass energy [GeV]" << endl;
        cout << "-minpt, -maxpt, -miny, -maxy" << endl;
        cout << "-ptstep stepsize" << endl;
        cout << "-pdf pdf" << endl;
        cout << "-fragfun ff" << endl;
        cout << "-nlo/-lo: use NLO/LO distributions" << endl;
        cout << "-x0 val: set x0 value for the BK solutions (overrides the value in BK file)" << endl;
        return 0;
        
    }


    double x0=-1;
    double y=-1;
    double sqrts=0;
    Order o;
    PDF* pdf;
    FragmentationFunction* ff;
    std::string datafile="";
    std::string datafile_probe="";
    Mode mode;
    Hadron final_particle;

    for (int i=1; i<argc; i++)
    {
        if (string(argv[i])=="-x0")
            x0 = StrToReal(argv[i+1]);
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
			ktfact_datafile2=argv[i+1];
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
        else if (string(argv[i]).substr(0,1)=="-")
        {
            cerr << "Unrecoginzed parameter " << argv[i] << endl;
            return -1;
        }
    }

    
    // Read data
    AmplitudeLib N(datafile);
    if (x0>0)
        N.SetX0(x0);
    if (y>=0)   // User set rapidity
        xbj=N.X0()*std::exp(-y);
    else // User set xbj
        y = std::log(N.X0()/xbj);

    N.InitializeInterpolation(xbj);
    cout << "# " << N.GetString() << endl;

    //**************** Different operation modes

    

    return 0;
}
    
    
