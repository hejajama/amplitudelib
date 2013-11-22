/* HERA f2 fitter
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012-2013
 * 
 * Computes F2 using the AmplitudeLib at the kinematical
 * values for which combined HERA F2 values are given
 * 
 * First argument is the BK solution datafile, 2nd argument is the 
 * HERA datafile with structure
 * Q^2 x inelasticity sigma_r error
 * Lines starting with # are ignored
 * inelasticity is the DIS y variable y=Q^2/(sx)
 * 
 * Prints the computed F2 values in the standard output, and the format is
 * Q^2 x y experimental_f2 f2_error theoryf2 (without normalization)
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
	if (argc != 3)
	{
		cout << "Syntax: " << argv[0] << " bksol  heradata " << endl;
		return 0;
	}
	string datafile = argv[1];
	string herafile = argv[2];
	
	cout << "# Reading BK equation solution from " << datafile << " and HERA data from " << herafile << endl;
	
	AmplitudeLib N(datafile);
	//N.SetX0(0.002);
	
	cout << "# " << N.GetString() << endl;
	cout <<"# x0 = " << N.X0() << endl;
	
	// Rear HERA data
	vector<double> xvals, yvals, qsqrvals, expvals, experrors;
	//NB: yvals is now vector of the inelasticities, y=Q^2/(sx)
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
		string x,qsqr,y,sigmar,err;
		stringstream l(line);
		l >> qsqr; l>>x; l>>y; l>>sigmar; l>>err;
		//if (std::abs(StrToReal(f2)) <1e-10) continue;
		if (StrToReal(x)>N.X0() or StrToReal(x)<1e-100) continue;
		qsqrvals.push_back(StrToReal(qsqr)); xvals.push_back(StrToReal(x));
		yvals.push_back(StrToReal(y));
		expvals.push_back(StrToReal(sigmar)); experrors.push_back(StrToReal(err)); 
	}
	file.close();

    VirtualPhoton wf;
    cout <<"# Wavefunction: " << wf << endl;
	
	cout << "# Q^2 [GeV^2]  x  y  HERA-\\sigma_r  HERA-err theory-\\sigma_r(light c b) " << endl;
	// Compute reduced cross section
	for (int i=0; i<xvals.size(); i++)
	{
		double x = xvals[i];
		const double mf=0.14;
		double xredef = x * (1.0 + 4.0*SQR(mf)/qsqrvals[i]);
		double sqrts = std::sqrt( qsqrvals[i]/(x * yvals[i]) );
		double sigmar_light = N.ReducedCrossSection(qsqrvals[i], std::log(N.X0()/x), sqrts, LIGHT);
		double sigmar_c = N.ReducedCrossSection(qsqrvals[i], std::log(N.X0()/x), sqrts, C);
		double sigmar_b = N.ReducedCrossSection(qsqrvals[i], std::log(N.X0()/x), sqrts, B);
		
		cout << qsqrvals[i] << " " << xvals[i] << " " << yvals[i] << " " << expvals[i] << " " << experrors[i] << " " << sigmar_light << " " << sigmar_c << " " << sigmar_b << endl;
	}
	
}
