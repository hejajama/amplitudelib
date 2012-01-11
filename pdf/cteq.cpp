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
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }
    if (q< 0.3 )
    {
        cerr << "q=" << q << " out of range at " << LINEINFO << endl;
        return 0;
    }
    // Codes
    int u=1; int d=2; int s=3; int c=4; int b=5; int g=0;
    int ubar=-1; int dbar=-2; 
    double result=0;
    // Antiquarks with minus sign 
    switch(p)
    {
        case U:
            result = ct10pdf_(u, x, q);
            break;
        case D:
            result = ct10pdf_(d,x,q);
            break;
        case UVAL:
            result = ct10pdf_(u,x,q) - ct10pdf_(ubar,x,q);
            break;
        case DVAL:
            result = ct10pdf_(d,x,q) - ct10pdf_(dbar,x,q);
            break;
        case USEA:
            result = ct10pdf_(ubar,x,q);
            break;
        case DSEA:
            result = ct10pdf_(dbar,x,q);
            break;
        case S:
            result = ct10pdf_(s,x,q);
            break;
        case C:
            result = ct10pdf_(c,x,q);
            break;
        case B:
            result = ct10pdf_(b,x,q);
            break;
        case G:
            result = ct10pdf_(g,x,q);
            break;
        default:
            cerr << "Parton " << p << " is not implemented " << LINEINFO << endl;            

    };

    return x*result;
}

// Default value of param is -1
void CTEQ::Initialize(int param)
{
    int set = param;
    if (param==-1)
        set=100;
    setct10_(set);
}

std::string CTEQ::GetString()
{
    return "CTEQ-TEA  CT10";
}