/*
 * Virtual class to hide different fragmentation functions
 * Uses directly dlib.f from http://www2.pv.infn.it/~radici/FFdatabase/
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

extern "C"
{
    REAL dlib_(REAL &z, REAL &Q2, REAL partons[11], int &ffset, int &fforder,
        int &ihadron, int &icharge, int &icp, int &ipi);
        // parton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar
}


// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
double FragmentationFunction::Evaluate(Parton p, Hadron h, REAL x, REAL q)
{
    cout <<" FragFun class is not tested, most likely doesn't work... " << LINEINFO << endl;
    return 0;
    /*
    REAL partons[11];

    int ihadron=0;
    int icharge=0;  // 0,1,2,3: 0 + - +&-
    int icp;    // 1: particle, 2: antiparticle, 3: sum (inactive for pi,h,p)
    
    if (h==HM)
    {
        ihadron = 3;
        icharge = 2;
    }
    else
    {
        cerr << "Parton/hadron combination is not supported... " << LINEINFO << endl;
        return 0;
    }

    int forder = 1;   // NLO
    int fset = 1;   // 1: K, 2: KKP, 3: BFGW
    int ipi = 1;    // best fit flag for BFGW

    REAL q2 = q*q;
    dlib_(x, q2, partons, fset, forder, ihadron, icharge, icp, ipi);

    if (p==U)
        return partons[6];
    else if (p==D)
        return partons[7];
    else if (p==G)
        return partons[5];
    else
    {
        cerr << "Parton " << p << " is not supported!  " << LINEINFO << endl;
        return 0;
    }
    */

}

FragmentationFunction::FragmentationFunction()
{
    order = NLO;

}

std::string FragmentationFunction::GetString()
{
    return "not specified";
}

void FragmentationFunction::SetOrder(Order o)
{
    order = o;
}

Order FragmentationFunction::GetOrder()
{
    return order;
}
