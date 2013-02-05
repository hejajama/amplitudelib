/*
 * PDF from UGD
 * Computes UGD from dipole amplitude and integrates it over k_T
 * to get PDF
 * 
 * Note: requires pointer to AmplitudeLib* and the scaling
 * factor \sigma_0/2
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#include "ugdpdf.hpp"
#include <string>

UGDPDF::UGDPDF()
{
	cerr << "UGDPDF requires AmplitdeLib*, panicking!" << endl;
	N=NULL;
	exit(1);
}

UGDPDF::UGDPDF(AmplitudeLib* N_, double sigma02_)
{
	N=N_; sigma02=sigma02_;
}

std::string UGDPDF::GetString()
{
	return "UGD PDF" ;
}

double UGDPDF::xq(double x, double q, Parton p)
{
	if (p != G)
	{
		//cerr << "UGDPDF can only compute gluon distribution! " << endl;
		return 0;
	}
	
	if (x < MinX())
	{
		cerr << "x < MinX()=" << MinX() << ", exit..." << endl;
		exit(1);
	}
	if (x>MaxX())
	{
		cerr << "x=" << x <<" > MaxX()=" << MaxX()<< endl;
		//x = MaxX();
		return 0;
		//exit(1);
	}
	
	return sigma02*N->xg(x, q);
}

UGDPDF::~UGDPDF()
{
	if (N != NULL)
		delete N;
}

double UGDPDF::MaxX()
{
	return N->X0();
}

double UGDPDF::MinX()
{
	return N->X0()*std::exp(-N->MaxY());
}

double UGDPDF::MinQ()
{
	return 0.1;
}
double UGDPDF::MaxQ()
{
	return 100000;
}
