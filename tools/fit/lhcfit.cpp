/* LHC hadronprod fitter
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012-2013
 * 
 * Computes hadron production pp collisions using AmplitudeLib 
 * 
 * First argument is the BK solution datafile, 2nd argument is the 
 * pp datafile with structure
 * miny maxy minpt maxpt hadron yield error
 * If minpt=maxpt (or y) then compute at fixed y/pt
 * Lines starting with # are ignored
 * 
 * Prints the computed F2 values in the standard output, and the format is
 * miny maxy minpt maxpt hadron sqrts yield error theory-yield (without normalization)
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
	cout << "# Hadron yield fitter" << endl;
	string datafile = argv[1];
	string lhcfile = argv[2];
	
	cout << "# Reading BK equation solution from " << datafile << " and hadron data from " << lhcfile << endl;
	
	AmplitudeLib N(datafile);
	AmplitudeLib N2(datafile);
	double s02 = 84.474/2.0;
	N.SetSigma02(s02);
	N2.SetSigma02(s02);
	
	// Reardata
	vector<double> minyvals,maxyvals,minptvals,maxptvals,yieldvals,errors,sqrtsvals; vector<Hadron> hadrons;
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
		string miny,maxy,minpt,maxpt,hadron,yield,err,sqrts;
		stringstream l(line);
		l >> miny; l>>maxy; l>>minpt; l>>maxpt; l>>hadron; l>>sqrts; l>>yield; l>>err;
		Hadron h;
		if (hadron=="CH") h=H;
		else if (hadron=="PI0") h=PI0;
		else if (hadron=="HM") h=HM;
		else { cerr << "Unknown hadron " << hadron << endl; exit(1); }
		
		minyvals.push_back(StrToReal(miny)); maxyvals.push_back(StrToReal(maxy));
		minptvals.push_back(StrToReal(minpt)); maxptvals.push_back(StrToReal(maxpt));
		sqrtsvals.push_back(StrToReal(sqrts));
		yieldvals.push_back(StrToReal(yield)); errors.push_back(StrToReal(err));
		hadrons.push_back(h);

	}
	file.close();
	
	cout <<"# sigma0 = " << N.Sigma02() << " GeV^-2" << endl;
	cout << "#  miny maxy minpt maxpt hadron sqrts yield error theory-yield   [1/GeV^2] " << endl;
	DSS ff;
	ff.SetOrder(LO);
	cout <<"# FF: " << ff.GetString() << endl;

	for (int i=0; i<minyvals.size(); i++)
	{
		double res = N.dHadronMultiplicity_dyd2pt_ktfact(minyvals[i], minptvals[i], sqrtsvals[i], &ff, hadrons[i], &N2);
		cout << minyvals[i] << " " << maxyvals[i] << " " << minptvals[i] << " " << maxptvals[i] << " " << hadrons[i] << "  " << sqrtsvals[i] << " " << yieldvals[i] << " " << errors[i] << " " << res << endl;
	}
	
}

