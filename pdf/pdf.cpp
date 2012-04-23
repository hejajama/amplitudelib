/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "pdf.hpp"

PDF::PDF()
{
	initialized=false;
}

void PDF::Initialize(int param)
{
}

std::string PDF::GetString()
{
    return "not specified";
}

PDF::~PDF()
{
}

void PDF::PlotPdf(double Q)
{
    cout << "# PDF at Q=" << Q << "GeV" << endl;
    cout <<"# x     up    down   s  gluon" << endl;
    for (double x=1e-4; x<1; x*=1.1)
    {
        cout << x << " " << xq(x,Q,U) << " " <<  xq(x,Q,D) << " "
        << xq(x,Q,S) << " " << xq(x,Q,G) << endl;
    }

}

double PDF::MaxX()
{
    cerr << "PDF::MaxX() is not implemented! " << LINEINFO << endl;
    return 1;
}
double PDF::MaxQ()
{
    cerr << "PDF::MaxQ() is not implemented! " << LINEINFO << endl;
    return 1e9;
}

double PDF::MinX()
{
    cerr << "PDF::MinX() is not implemented! " << LINEINFO << endl;
    return 0;
}
double PDF::MinQ()
{
    cerr << "PDF::MinQ() is not implemented! " << LINEINFO << endl;
    return 0;
}
