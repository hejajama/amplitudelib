/*
 * Virtual class to hide different fragmentation functions
 * Uses directly dlib.f from http://www2.pv.infn.it/~radici/FFdatabase/
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
 */

#include <string>
#include <sstream>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

extern "C"
{
    double dlib_(double &z, double &Q2, double partons[11], int &ffset, int &fforder,
        int &ihadron, int &icharge, int &icp, int &ipi);
        // parton = 5,4,3,2,1,0,-1,...,-5 means b,c,s,d,u,g,ubar,...,bbar
}

using namespace Amplitude;

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
double FragmentationFunction::Evaluate(Parton p, Hadron h, double x, double q)
{
    cout <<" FragFun class is not tested, most likely doesn't work... " << LINEINFO << endl;
    return 0;
    /*
    double partons[11];

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

    double q2 = q*q;
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

void FragmentationFunction::Test()
{
	cerr << "FragmentationFunction::Test() is not implemented" << endl;
}


std::string ParticleStr(Hadron h)
{
	switch(h)
	{
		case P:
			return "proton";
		case PI0:
			return "pi0";
		case H:
			return "charged hadron";
		case HM:
			return "h-";
		case HP:
			return "h+";
		case PIP:
			return "pi+";
		case PIM:
			return "pi-";
		default:
			std::stringstream ss;
			ss << "unkown particle type " << h;
			return ss.str();
	}
	
	return "ERROR!";
	
}
