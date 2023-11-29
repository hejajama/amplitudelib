/*
 * DSS Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */
#include <cassert>
#include "dss.hpp"
#include "../tools/config.hpp"
#include "fragmentation.hpp"
#include <cstdlib>
#include "../tools/tools.hpp"


extern "C"
{
    // ih: input hadron, 1=pion, 2=kaon, 3=proton, 4=charged hadrons
    // ic: charge, 0=0, 1=+, -1=-
    // parton: 
		//0    1    2    3    4    5    6    7    8     9    10
		//g    u   ubar  d   dbar  s   sbar  c   cbar   b   bbar
    // io: order, 0=LO, 1=NLO
    // NOTE: cbar and bbar not in DSS, returns c or b instead
	// result is the fragmentation function (X factor divided out)
    
    // version: 1 = DSS07, 2=DSS14 (only pions)
    void fdss_ (int &version,int &ih, int &ic, int &io, double &x, double& q2, double& u,
        double &ub, double &d, double &db, double &s, double &sb, double &c,
        double &b, double &g);
    /*
    void fdss_(int& hadron, int& charge, int& order, double& z,
              double& scalesqr, int& parton, double& result);
              */

    extern struct{
        double fini;
    } fragini_;
}


using namespace Amplitude;

// D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
double DSS::Evaluate(Parton p, Hadron h, double x, double qs)
{
    if (!(h == PI0 or h == PIP or h==PIM) and version==DSS14)
    {
        std::cerr << "DSS14 does not support anything else than pions" << std::endl;
        exit(1);
    }

    if (x<0.05 or x>1)
    {
        cerr << "z=" << x <<" out of range [0.05, 1] at " << LINEINFO << endl;
        return 0;
    }
    if (qs*qs<1 or qs*qs>1e5)
    {
		cerr << "Q=" << qs << " GeV out of Q^2 range [1,1e5] GeV^2. " << LINEINFO << endl;
		if (qs*qs<1)
			qs=1;
		else
			qs=std::sqrt(1e5);
        
        //exit(1);
        //return 0;
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
		ih=4; ic=-1;
    } else if (h==HP)
    {
		ih=4; ic=1;
    } 
    else if (h==H) // charged hadron
    {
		return Evaluate(p, HP, x, qs) + Evaluate(p, HM, x, qs);
	}
    else if (h==PI)	// pi^+ + pi^-
    {    
		return Evaluate(p, PIP, x, qs) + Evaluate(p, PIM, x, qs);
	}
    else
    {
        cerr << "Hadron " << h << "(" << ParticleStr(h) <<") is not supported by DSS " << LINEINFO << endl;
        return 0;
    }

    if (!initialized)
        fragini_.fini = 0;
    initialized=true;

    double qsqr = qs*qs;
    //double u, ubar, d, dbar, s, sbar, c, b, g;
    //fdss_(ih, ic, io, x, qsqr, u, ubar, d, dbar, s, sbar, c, b, g);
    //double result;
    // New DSS14 wrapper partons: U, UB, D, DB, S, SB,   C,           B,       GL     
    int parton=0;
    switch(p)
    {
        case U:
            //result = u/x;
            parton=0;
            break;
        case D:
			parton=2;
            //result = d/x;
            break;
        case UBAR:
			parton=1;
            //result = ubar/x;
            break;
        case DBAR:
			parton=3;
            //result = dbar/x;
            break;
        case S:
			parton=4;
            //result = s/x;
            break;
        case SBAR:
            parton=5;
            //result = sbar/x;
            break;
        case C:
            parton = 6;
            break;
		case CBAR:
			//parton = 8;
            parton=6;
			break;
        case B:
            //parton = 9;
            parton=7;
            break;
		case BBAR:
			//parton=10;
            parton=7;
			break;
        case G:
			parton=8;
            break;
        default:
            cerr << "Parton " << PartonToString(p) << " is not supported! " << LINEINFO << endl;
            return 0;
    }

    double output[9];
    int v=-1;
    if (version==DSS07) v=1;
    else if (version==DSS14) v=2;
    else
    {
        cerr << "Uknown DSS version! " << endl;
        exit(1);
    }
    fdss_(v, ih, ic, io, x, qsqr, output[0],output[1],output[2],output[3],output[4],output[5],
        output[6],output[7],output[8] );
	//fdss_(ih, ic, io, x, qsqr, parton, result);
    
    return output[parton]/x; // this returns D(z), not z*D(z) 
}



std::string DSS::GetString()
{
    if (order==LO) return "DSS LO";
    return "DSS NLO";
}

DSS::DSS()
{
    initialized = false;
    fragini_.fini=0;
    version=DSS07;
}
 
void DSS::SetOrder(Order o)
{
	order=o;
	initialized=false;
	fragini_.fini=0;
}

/*
 * Run tests
 * "Correct" results are obtained using FFgenerator
 *  http://lapth.in2p3.fr/ffgenerator/
 */
void DSS::Test()
{
	
	cout <<"#----- Testing DSS FF" << endl;
	cout << "#NLO:" << endl;
	Order orig_order = order;
	
    double accuracy=1e-2;
    
	SetOrder(NLO);
	double result, cor;
	result=0.1*Evaluate(U, PI0, 0.1, std::sqrt(10)); cor=0.64450;
	cout <<"z*D_{u->pi0}(Q^2=10GeV^2, z=0.1) = " << result << " (correct " << cor << ")" <<  endl;
    assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.2*Evaluate(D, PI, 0.2, std::sqrt(50)); cor=0.97394;
	cout <<"z*D_{d->pi+ pi-}(Q^2=50GeV^2, z=0.2) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.3*Evaluate(S, HM, 0.3, std::sqrt(20)); cor=0.64972;
	cout <<"z*D_{s->h-}(Q^2=20GeV^2, z=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.3*Evaluate(G, PI0, 0.3, std::sqrt(20)); cor=0.37437;
	cout <<"z*D_{g->pi0}(Q^2=20GeV^2, z=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result = 0.05 * Evaluate(U, H, 0.05, std::sqrt(20)); cor=1.7988;
	cout << "z*D_{u->charged hadron}(Q^2=20 GeV^2, z=0.05) = " << result << " (correct " << cor << ")" << endl;
	assert(std::abs(result-cor)/cor<accuracy);
		
	cout <<"#LO:" << endl;
	SetOrder(LO);
	
		result=0.1*Evaluate(U, PI0, 0.1, std::sqrt(10)); cor=0.64300;
	cout <<"z*D_{u->pi0}(Q^2=10GeV^2, x=0.1) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.2*Evaluate(D, PI, 0.2, std::sqrt(50)); cor=1.0146;
	cout <<"z*D_{d->pi+ pi-}(Q^2=50GeV^2, x=0.2) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.3*Evaluate(S, HM, 0.3, std::sqrt(20)); cor=0.74499;
	cout <<"z*D_{s->h-}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result=0.3*Evaluate(G, PI0, 0.3, std::sqrt(20)); cor=0.6663;
	cout <<"z*D_{g->pi0}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	result = 0.05 * Evaluate(U, H, 0.05, std::sqrt(20)); cor=1.9169;
	cout << "z*D_{u->charged hadron}(Q^2=20 GeV^2, z=0.05) = " << result << " (correct " << cor << ")" << endl;
	assert(std::abs(result-cor)/cor<accuracy);
	
	
	cout << "All tests done, if no errors were shown, all tests passed!" << endl;
	SetOrder(orig_order);
	
}
