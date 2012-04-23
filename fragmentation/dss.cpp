/*
 * DSS Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "dss.hpp"
#include "../tools/config.hpp"
#include "fragmentation.hpp"

extern "C"
{
    // ih: input hadron, 1=pion, 2=kaon
    // ic: charge, 0=0, 1=+, -1=-
    // io: order, 0=LO, 1=NLO
    
    void fdss_ (int &ih, int &ic, int &io, double &x, double& q2, double& u,
        double &ub, double &d, double &db, double &s, double &sb, double &c,
        double &b, double &g);

    extern struct{
        int fini;
    } fragini_;
}

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
double DSS::Evaluate(Parton p, Hadron h, double x, double qs)
{
    if (x<0.05 or x>1)
    {
        cerr << "x=" << x <<" out of range [0.05, 1] at " << LINEINFO << endl;
        return 0;
    }
    if (qs*qs<1 or qs*qs>1e5)
    {
        cerr << "Q=" << qs << " GeV out of Q^2 range [1,1e5] GeV^2. " << LINEINFO << endl;
        return 0;
    }

    int ih, ic, io;
    if (order == NLO) io=1;
    else io=0;

    if (h == PIM)
    {
        ih=1; ic=-1;
    } else if (h==PIP)
    {
        ih=1; ic=1;
    } else if (h==PI0)
    {
        ih=1; ic=0;
    } else if (h==KP)
    {
        ih=2; ic=1;
    } else if (h==KM)
    {
        ih=2; ic=-1;
    } else if (h==K0)
    {
        ih=2; ic=0;
    }
    else if (h==HM)
    {
        return Evaluate(p, PIM, x, qs) + Evaluate(p, KM, x, qs);
    } else if (h==HP)
    {
        return Evaluate(p, PIP, x, qs) + Evaluate(p, KP, x, qs);
    }     
    else
    {
        cerr << "Hadron " << h << " is not supported by DSS " << LINEINFO << endl;
        return 0;
    }

    if (!initialized)
        fragini_.fini = 0;
    initialized=true;

    double qsqr = qs*qs;
    double u, ubar, d, dbar, s, sbar, c, b, g;
    fdss_(ih, ic, io, x, qsqr, u, ubar, d, dbar, s, sbar, c, b, g);
    double result;
    switch(p)
    {
        case U:
            result = u/x;
            break;
        case D:
            result = d/x;
            break;
        case UBAR:
            result = ubar/x;
            break;
        case DBAR:
            result = dbar/x;
            break;
        case S:
            result = s/x;
            break;
        case SBAR:
            result = sbar/x;
            break;
        case G:
            result = g/x;
            break;
        default:
            cerr << "Parton " << p << " is not supported! " << LINEINFO << endl;
            return 0;
    }

    
    return result;  // Shoudn't end up here
}



std::string DSS::GetString()
{
    return "DSS";
}

DSS::DSS()
{
    initialized = false;
    fragini_.fini=0;
}
 
