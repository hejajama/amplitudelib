/* LHC hadronprod fitter
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 * 
 * Computes hadron production in lhc kinematics using the AmplitudeLib 
 * 
 * First argument is the BK solution datafile, 2nd argument is the 
 * LHC datafile with structure
 * y hadron pt yield error
 * Lines starting with # are ignored
 * 
 * Prints the computed F2 values in the standard output, and the format is
 * y hadron sqrts pt yield error theory-yield (without normalization)
 */
 


#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include <amplitudelib/datafile.hpp>
#include <fragmentation/fragmentation.hpp>
#include <fragmentation/dss.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <gsl/gsl_errno.h>
using namespace std;

int main(int argc, char* argv[])
{
	gsl_set_error_handler(&ErrHandler);
	cout << "# LHC fitter" << endl;
	string datafile = argv[1];
	string lhcfile = argv[2];
	
	cout << "# Reading BK equation solution from " << datafile << " and LHC data from " << lhcfile << endl;
	
	AmplitudeLib N(datafile);
	AmplitudeLib N2(datafile);
	
	// Rear HERA data
	vector<double> yvals,ptvals,yieldvals,errors,sqrtsvals; vector<Hadron> hadrons;
    ifstream file(lhcfile.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << lhcfile.c_str() << endl;
        return -1;
    }
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#" or line.length()<10)
			continue;
		string y,hadron,pt,yield,err,sqrts;
		stringstream l(line);
		l >> y; l>>hadron; l>>sqrts; l>>pt; l>>yield; l>>err;
		Hadron h;
		if (hadron=="CH") h=H;
		else if (hadron=="PI0") h=PI0;
		else { cerr << "Unknown hadron " << hadron << endl; exit(1); }
		
		yvals.push_back(StrToReal(y)); ptvals.push_back(StrToReal(pt));
		sqrtsvals.push_back(StrToReal(sqrts));
		yieldvals.push_back(StrToReal(yield)); errors.push_back(StrToReal(err));
		hadrons.push_back(h);

	}
	file.close();
	
	cout << "# y hadron  sqrts pt yield err theory  [1/GeV^2] " << endl;
	DSS ff;
	ff.SetOrder(LO);
	// Compute F2
	for (int i=0; i<yvals.size(); i++)
	{
		double res = N.dHadronMultiplicity_dyd2pt_ktfact(yvals[i], ptvals[i], sqrtsvals[i], &ff, hadrons[i], &N2);
		cout << yvals[i] << " " << hadrons[i] << "  " << sqrtsvals[i] << " " << ptvals[i] << " " << yieldvals[i] << " " << errors[i] << " " << res << endl;
	}
	
}

