/*
 * Wrapper calss for CTEQ-TEA parton distribution functions CT10
 * Real code is in file CT10Pdf.f downloaded from
 * http://hep.pa.msu.edu/cteq/public/index.html
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "cteq.hpp"
#include <string>

double CTEQ::xq(double x, double q, Parton p)
{
	if (!initialized)
	{
		cerr <<" CTEQ parton distribution function is not initialized, "
			<< "can't evaluate CTEQ::xq() at " << LINEINFO << endl;
		return 0;
	}
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }
    double minq=0.3;
    if (order==LO) minq=1.3;
    if (q< minq )
    {
        cerr << "q=" << q << " out of range at " << LINEINFO << endl;
        return 0;
    }
    // Codes
    int u=1; int d=2; int s=3; int c=4; int b=5; int g=0;
    int ubar=-1; int dbar=-2; 
    double result=0;
    // Antiquarks with minus sign 
    
    // LO or NLO
    double (*f)(int& iparton, double& x, double& q);
    if (order==LO) f = ctq6pdf_;
    else if (order==NLO) f = ct10pdf_;
    
    switch(p)
    {
        case U:
            result = f(u, x, q);
            break;
        case D:
            result = f(d,x,q);
            break;
        case UVAL:
            result = f(u,x,q) - f(ubar,x,q);
            break;
        case DVAL:
            result = f(d,x,q) - f(dbar,x,q);
            break;
        case USEA:
            result = f(ubar,x,q);
            break;
        case DSEA:
            result = f(dbar,x,q);
            break;
        case S:
            result = f(s,x,q);
            break;
        case C:
            result = f(c,x,q);
            break;
        case B:
            result = f(b,x,q);
            break;
        case G:
            result = f(g,x,q);
            break;
        default:
            cerr << "Parton " << p << " is not implemented " << LINEINFO << endl;            

    };

    return x*result;
}

// Default value of param is -1
void CTEQ::Initialize(int param)
{
	// Useless...
}

void CTEQ::SetOrder(Order o)
{
	order=o;
	
	if (order==NLO)
	{
		int set=100;
		setct10_(set);
	}
	else if (order==LO)
	{
		int set = 3;
		setctq6_(set);
	}
	initialized=true;

}	

std::string CTEQ::GetString()
{
	if (order==LO)
		return "CTEQ6 LO";
	
    return "CTEQ-TEA NLO  CT10";
}

double CTEQ::MinQ()
{
    return 0.3;
}

double CTEQ::MaxQ()
{
    return 100000;
}

CTEQ::CTEQ()
{
	SetOrder(NLO);		// Default
}
