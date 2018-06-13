/*
 * Wrapper calss for EPS09 nuclear parton distribution functions 
 * Real code is in file EPS09.f downloaded from
 * https://www.jyu.fi/fysiikka/en/research/highenergy/urhic/eps09
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#include "eps09.hpp"
#include <string>
#include <cstdlib>

using namespace Amplitude;

double EPS09::xq(double x, double q, Parton p)
{
	if (!initialized)
	{
		cerr <<" EPS parton distribution function is not initialized, "
			<< "can't evaluate EPS09::xq() at " << LINEINFO << endl;
		return 0;
	}
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }
    if (q< MinQ() or q>MaxQ())
    {
		cerr << "q=" << q << " out of range at " << LINEINFO << endl;
		return 0;
    }
    // Codes
    double ruv, rdv, ru, rd, rs, rc, rb, rg;
    double modification=0;
    int set=1;
    int o;
    if (order==LO)
		o=1;
	else
		o=2;
    
    double Z=-1;
    if (A == 208)
        Z=82;
    else
    {
        cerr << "Unknown nucleus " << A << "! " << LINEINFO << endl;
        return 0;
    }
	
	eps09_(o, set, A, x, q, ruv, rdv, ru, rd, rs, rc, rb, rg);
    
    switch(p)
    {
        case U:
            return A*(Z/A * ru*cteq.xq(x,q,U) + (A-Z)/A * rd*cteq.xq(x,q,D));
        case UBAR:
            return A*(Z/A * ru*cteq.xq(x,q,UBAR) + (A-Z)/A * rd*cteq.xq(x,q,DBAR));
        case D:
            return A*(Z/A * rd*cteq.xq(x,q,D) + (A-Z)/A * ru*cteq.xq(x,q,U));
        case DBAR:
            return A*(Z/A * rd*cteq.xq(x,q,DBAR) + (A-Z)/A * ru*cteq.xq(x,q,UBAR));
        case S:
        case SBAR:
            return A*rs*cteq.xq(x,q,p);
        case C:
        case CBAR:
            return A*rc*cteq.xq(x,q,p);
        case B:
        case BBAR:
            return A*rb*cteq.xq(x,q,p);
        case G:
            return A*cteq.xq(x, q, p) * rg;
        default:
            cerr << "Parton " << p << " is not implemented " << LINEINFO << endl;            

    };

    return 0;
    
}

// Default value of param is -1
void EPS09::Initialize(int param)
{
	cteq.Initialize();
}

void EPS09::SetOrder(Order o)
{
	order=o;
	cteq.SetOrder(o);
	Initialize();
	initialized=true;

}	

std::string EPS09::GetString()
{
	if (order==LO)
		return "EPS09 LO";
    return "EPS09 NLO";
}

double EPS09::MinQ()
{
	return 1.3;	
}

double EPS09::MaxQ()
{
    return 1000;
}

EPS09::EPS09()
{
	SetOrder(NLO);		// Default
}

/*
 * Run tests,
 * Check that we obtain correct numbers in some special cases
 */
void EPS09::Test()
{
	
}

int EPS09::SetA(int A_)
{
	A=A_;
	return 0;
}
