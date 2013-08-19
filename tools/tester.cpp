/*
 * Tester program for AmplitudeLib
 */

#include "config.hpp"
#include "tools.hpp"
#include "../amplitudelib/amplitudelib.hpp"
#include "../amplitudelib/datafile.hpp"
#include "../pdf/pdf.hpp"
#include "../pdf/cteq.hpp"
#include "../fragmentation/dss.hpp"
#include <gsl/gsl_errno.h>
#include <string>

using namespace std;

int main()

{
	gsl_set_error_handler(&ErrHandler);
	cout << "This is a tester program for AmplitudeLib" << endl;
	cout << "=========================================" << endl;
	
	string dataf = "mv1_noevolution_qs02.dat";
	cout << "Loading data from " << dataf << endl << "assuming that at all "
		<<"rapidities it is just MV-model with Qs0^2=0.2GeV^2" << endl;
		
	AmplitudeLib N(dataf, false);
	double y=0.1;
	double res,correct;
	N.InitializeInterpolation(y);
	
	cout << "===== TEST: evaluation of N(r) ===== " << endl;
	correct=0.091845; res=N.N(1,y);
	cout << "N(r=1) = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct; else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: 2d FT of S(r)=1-N(r) ===== " << endl;
	correct=1.17962; res=N.S_k(1.0, y);
	cout << "FT S(k=1) = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: 2d FT of S^2(r)=(1-N(r))^2 ===== " << endl;
	correct=3.31528; res=N.S_k(1, y, false, 2.0);
	cout << "FT S(k=1)^2 = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: 2d FT of N(r)/2\\pi r^2 ===== " << endl;
	correct=0.00958729; res=N.N_k(2, y);
	cout << "FT N(k=2)/(2pi r^2) = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: d ln N(r) / d ln r ===== " << endl;
	correct=0.733868; res=N.LogLogDerivative(2, y);
	cout << "d ln N(r=2)/d ln r = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST UGD  ===== " << endl;
	N.SetRunningCoupling(FIXED);
	correct=0.0891024; res=N.UGD(1, y, 1, 1);
	cout << "UGD(x=x0, Q^2=1 GeV^2)=" << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST UGD -> xg ===== " << endl;
	N.SetRunningCoupling(FIXED);
	correct=0.266333; res=N.xg(N.X0(), 4);
	cout << "xg(x=x0, Q^2=16 GeV^2)=" << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	
	
	cout << "===== TESTING INTERPOLATOR =====" << endl;
	{
	std::vector<double> y; y.push_back(0); y.push_back(1); y.push_back(4);  y.push_back(9); y.push_back(16);
	std::vector<double> x; x.push_back(0); x.push_back(1); x.push_back(2);  x.push_back(3); x.push_back(4);
	Interpolator interp(x,y); interp.Initialize();
	correct=6.25; res = interp.Evaluate(2.5);	// 2.5^2=6.25
	cout << "6.25^2 = " << res; if (abs(res-correct)/correct>0.01) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	// derivative at x=1.5 = 2*1.5=3
	correct=3; res=interp.Derivative(1.5);
	cout << "D(x^2, x=1.5) = " << res; if (abs(res-correct)/correct>0.05) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	}
	
	
	cout << "===== TEST: CTEQ pdf ===== " << endl;
	CTEQ pdf;
	pdf.Initialize();
	pdf.Test();
	
	cout << "===== TEST: DSS FF ===== " << endl;
	DSS ff;
	ff.Test();
	
		
		
		
		
	return 0;
	
	
}
