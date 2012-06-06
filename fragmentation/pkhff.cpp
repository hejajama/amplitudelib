/*
 * PKHFF Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "pkhff.hpp"
#include "../tools/config.hpp"
#include "fragmentation.hpp"

extern "C"
{
    // parton labels in output array:
    // output arrays are u[1] = u-quark, u[2]=ubar-quark
    
    void pkhff_(int& set, int& charge, REAL& x, REAL &qsqr,
        REAL u[2], REAL d[2], REAL s[2], REAL c[2], REAL b[2], REAL& g, int &ifini);
}

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
REAL PKHFF::Evaluate(Parton p, Hadron h, REAL x, REAL qs)
{
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }

    int charge=0; int set=0;
    if (h == HP)
    {
        charge=1;
        set=6;
    } else if (h==HM)
    {
        charge=2;
        set=6;
    } else if (h==H)
    {
        charge=3;
        set=6;
    } else if (h==PIP)
    {
        charge=1;
        set=2;
    }
    else if (h==PIM)
    {
        charge=2;
        set=2;
    }
    else if (h==PI0)
    {
        return 0.5*(Evaluate(p, PIP, x, qs) + Evaluate(p, PIM, x, qs));
    }
     else if (h==PI)
    {
        charge=3;
        set=2;
    }
    else
    {
        cerr << "Hadron " << h << " is not supported by the PKHFF fragmentation! "
            << LINEINFO << endl;
        return 0;
    }
        

    REAL u[2], d[2], s[2], c[2], b[2];
    REAL g=0;
    REAL qsqr=qs*qs;
    if (qsqr < 1.01 or qsqr>5.99e6)
    {
        cerr << "Q^2 out of range for PKHFF fragmentation! " << LINEINFO << endl;
        return 0;
    }
    pkhff_(set, charge, x, qsqr, u, d, s, c, b, g, ifini);

    REAL result=0;
    if (p==U)
        result=u[0];
    else if (p==D)
        result=d[0];
    else if (p==S)
        result = s[0];
    else if (p==C)
        result = c[0];
    else if (p==B)
        result = b[0];
    else if (p==G)
        result = g;
    else
    {
        cerr << "Parton " << p<< " is not supported by PKHFF fragmentation! "
            << LINEINFO << endl;
        return 0;
    }
        
    
    return result;

}



std::string PKHFF::GetString()
{
    return "PKHFF";
}

PKHFF::PKHFF()
{
    ifini=0;
    cerr << "PHKFF is not well tested/may not support NLO/LO" << endl;

}
 
