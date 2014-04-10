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
#include "interpolation2d.hpp"
#include <gsl/gsl_errno.h>
#include <string>

using namespace std;
using namespace Amplitude;

int main()

{
	gsl_set_error_handler(&ErrHandler);
	cout << "This is a tester program for AmplitudeLib" << endl;
	cout << "=========================================" << endl;
	
	string dataf = "mv1_noevolution_qs02.dat";
	cout << "Loading data from " << dataf << endl << "assuming that at all "
		<<"rapidities it is just MV-model with Qs0^2=0.2GeV^2" << endl;
		
	AmplitudeLib N(dataf, false);
    N.GetQCD().SetRunningCoupling(FIXED);
	double y=0.1;
    double xbj = N.X0()*std::exp(-y);
	double res,correct;
	N.InitializeInterpolation(xbj);
	
	cout << "===== TEST: evaluation of N(r) ===== " << endl;
	correct=0.091845; res=N.N(1,xbj);
	cout << "N(r=1) = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct; else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: 2d FT of S(r)=1-N(r) ===== " << endl;
	correct=1.17962; res=N.S_k(1.0, xbj);
	cout << "FT S(k=1) = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: 2d FT of S^2(r)=(1-N(r))^2 ===== " << endl;
	correct=3.31528; res=N.S_k(1, xbj, FUNDAMENTAL, 2.0);
	cout << "FT S(k=1)^2 = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST WW UGD =  2d FT of N(r)/2\\pi r^2 ===== " << endl;
	correct=0.00958729; res=N.WW_UGD(2, xbj);
	cout << "FT N(k=2)/(2pi r^2) [WW ugd] = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST: d ln N(r) / d ln r ===== " << endl;
	correct=0.733868; res=N.LogLogDerivative(2, xbj);
	cout << "d ln N(r=2)/d ln r = " << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST Dipole UGD  ===== " << endl;
	//N.SetRunningCoupling(FIXED);
	correct=0.0891024; res=N.Dipole_UGD(1, xbj);
	cout << "UGD(x=x0, Q^2=1 GeV^2)=" << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	cout << "===== TEST UGD -> xg ===== " << endl;
	correct=0.266333; res=N.xg(xbj, 4);
	cout << "xg(x=x0, Q^2=16 GeV^2)=" << res; if (abs(res-correct)/correct>0.001) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	
	
	
	cout << "===== TESTING INTERPOLATOR =====" << endl;
	{
	std::vector<double> y; y.push_back(0); y.push_back(1); y.push_back(4);  y.push_back(9); y.push_back(16);
	std::vector<double> x; x.push_back(0); x.push_back(1); x.push_back(2);  x.push_back(3); x.push_back(4);
	Interpolator interp(x,y); interp.Initialize();
	correct=6.25; res = interp.Evaluate(2.5);	// 2.5^2=6.25
	cout << "2.5^2 = " << res << " (correct: " << correct << ") "; if (abs(res-correct)/correct>0.01) cout << " TEST FAILED!" ;else cout << " OK!";
	cout << endl;
	
	// derivative at x=1.5 = 2*1.5=3
	correct=3; res=interp.Derivative(1.5);
	cout << "D(x^2, x=1.5) = " << res; if (abs(res-correct)/correct>0.05) cout << " TEST FAILED! Correct: " << correct;else cout << " OK!";
	cout << endl;
	}

    cout << "===== TESTING 2D INTERPOLATOR =====" << endl;
	{
	std::vector<double> axes; axes.push_back(0); axes.push_back(1); axes.push_back(2); axes.push_back(3); axes.push_back(4);
    // Test function x*y^2
    // Row 1. x=0
    std::vector<double> row1; row1.push_back(0); row1.push_back(0);  row1.push_back(0);  row1.push_back(0);  row1.push_back(0);
    // Row 2, x=1
    std::vector<double> row2; row2.push_back(0); row2.push_back(1);  row2.push_back(4);  row2.push_back(9);  row2.push_back(16);
    // Row 3, x=2
    std::vector<double> row3; row3.push_back(0); row3.push_back(2);  row3.push_back(8);  row3.push_back(18);  row3.push_back(32);
    // Row 4, x=3
    std::vector<double> row4; row4.push_back(0); row4.push_back(3);  row4.push_back(12);  row4.push_back(27);  row4.push_back(48);
    // Row 4, x=4
    std::vector<double> row5; row2.push_back(0); row5.push_back(4);  row5.push_back(16);  row5.push_back(36);  row5.push_back(64);
    std::vector<std::vector<double> > data;
    data.push_back(row1); data.push_back(row2); data.push_back(row3); data.push_back(row4); data.push_back(row5);
    
    
	Interpolator2D interp2(axes, data);
    double interpx=1.2; double interpy=2.8;
    correct = interpx*interpy*interpy;
	res = interp2.Evaluate(interpx, interpy);
	cout << "x*y^2 (x=" << interpx <<", y=" << interpy <<") = " << res << " (correct: " << correct << ")" ; if (abs(res-correct)/correct>0.05) cout << " TEST FAILED!";else cout << " OK!";
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
