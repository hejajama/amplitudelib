/* HERA f2 fitter
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 * 
 * Computes F2 using the AmplitudeLib at the kinematical
 * values for which combined HERA F2 values are given
 * 
 * First argument is the BK solution datafile, 2nd argument is the 
 * HERA datafile with structure
 * Q^2 x experimental_f2 f2_error
 * Lines starting with # are ignored
 * 
 * Prints the computed F2 values in the standard output, and the format is
 * Q^2 x experimental_f2 f2_error theoryf2 (without normalization)
 */
 


#include <tools/config.hpp>
#include <tools/tools.hpp>
#include <amplitudelib/amplitudelib.hpp>
#include <amplitudelib/datafile.hpp>
#include <string>
#include <fstream>
#include <sstream>
#include <gsl/gsl_errno.h>
using namespace std;

int main(int argc, char* argv[])
{
	gsl_set_error_handler(&ErrHandler);
	cout << "# F2 fitter" << endl;
	string datafile = argv[1];
	string herafile = argv[2];
	
	cout << "# Reading BK equation solution from " << datafile << " and HERA data from " << herafile << endl;
	
	AmplitudeLib N(datafile);
	
	// Rear HERA data
	vector<double> xvals, qsqrvals, expvals, experrors;
    ifstream file(herafile.c_str());
    if (!file.is_open())
    {
        cerr << "ERROR! Coudn't read file " << file << endl;
        return -1;
    }
    
    while(!file.eof() )
    {
        string line;
        getline(file, line);
        if (line.substr(0, 1)=="#")
			continue;
		string x,qsqr,f2,f2err;
		stringstream l(line);
		l >> qsqr; l>>x; l>>f2; l>>f2err;
		if (std::abs(StrToReal(f2)) <1e-10) continue;
		qsqrvals.push_back(StrToReal(qsqr)); xvals.push_back(StrToReal(x));
		expvals.push_back(StrToReal(f2)); experrors.push_back(StrToReal(f2err)); 
	}
	file.close();
	
	cout << "# Q^2 [GeV^2]  x  HERA-F2  HERA-err Theory-f2 (no normalization) " << endl;
	// Compute F2
	for (int i=0; i<xvals.size(); i++)
	{
		double x = log(N.X0()/xvals[i]);
		double f2 = N.F2(qsqrvals[i], x);
		cout << qsqrvals[i] << " " << xvals[i] << " " << expvals[i] << " " << experrors[i] << " " << f2 << endl;
	}
	
}
