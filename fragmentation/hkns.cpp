/*
 * KKP Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "hkns.hpp"
#include "../tools/config.hpp"
#include "fragmentation.hpp"

extern "C"
{
 /*     FF(I) --> HKNS fragmentation functions (FFs).
       I = -5 ... b-bar quark (D_b-bar = D_b)
           -4 ... c-bar quark (D_c-bar = D_c)
           -3 ... s-bar quark
           -2 ... d-bar quark 
           -1 ... u-bar quark 
            0 ... gluon D_g(x)
            1 ... u quark 
            2 ... d quark 
            3 ... s quark 
            4 ... c quark 
            5 ... b quark
            */
    
    void hknsff_(double &q2, double &x, int &iset, int &icharge, double partons[11], double grad[17][11]);
}

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
REAL HKNS::Evaluate(Parton p, Hadron h, REAL x, REAL qs)
{
    if (x<0 or x>1)
    {
        cerr << "x=" << x <<" out of range at " << LINEINFO << endl;
        return 0;
    }
    if (qs<1 or qs>1e4)
    {
        cerr << "q=" << qs << " GeV out of range [1,1e4] GeV at " << LINEINFO << endl;
        return 0;
    }

    
    int charge=0; int set=0;

    if (h == PIP)
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
        return 0.5*( Evaluate(p, PIP, x, qs) + Evaluate(p, PIM, x, qs) );
    }
    else if (h==PI)
    {
        charge=3;
        set=2;
    }
    else if (h==PP)
    {
        charge=1;
        set=6;
    }
    else if (h==PM)
    {
        if (p==U)
            return Evaluate(UBAR, PP, x, qs);
        else if (p==D)
            return Evaluate(DBAR, PP, x, qs);
        else if (p==S)
            return Evaluate(SBAR, PP, x, qs);
        else if (p==G)
            return Evaluate(G, PP, x, qs);
        else
        {
            cerr << "Can't calculate p^- with given quark... " << LINEINFO << endl;
            return 0;
        }
    }
    else if (h==KP)
    {
        charge=1;
        set=4;
    }
    else if (h==KM)
    {
        charge=2;
        set=4;
    }
    else if (h==HM)
    {
        return Evaluate(p, PIM, x, qs) + Evaluate(p, KM, x, qs) + Evaluate(p, PM, x, qs);
    }
    else
    {
        cerr << "Hadron " << h << " is not supported by the HKNS fragmentation! "
            << LINEINFO << endl;
        return 0;
    }
    double qs2 = qs*qs;
    double errors[17][11];
    REAL partons[11];
    hknsff_(qs2, x, set, charge, partons, errors);
    REAL result=0;
    
    if (p==G) result = partons[5];
    else if (p==U) result = partons[6];
    else if (p==D) result = partons[7];
    else if (p==S) result = partons[8];
    else if (p==UBAR) result = partons[4];
    else if (p==DBAR) result = partons[3];
    else if (p==SBAR) result = partons[2];
    else
        cerr << "Parton " << p << " is not supported! " << LINEINFO << endl;
    
    return result;

}



std::string HKNS::GetString()
{
    return "HKNS";
}

HKNS::HKNS()
{
	cerr << "HKNS is not well tested and may not support LO/NLO!" << endl;
}
 
