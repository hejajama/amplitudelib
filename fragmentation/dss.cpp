/*
 * DSS Fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include "dss.hpp"
#include "../tools/config.hpp"
#include "fragmentation.hpp"
#include <cstdlib>

extern "C"
{
    // ih: input hadron, 1=pion, 2=kaon
    // ic: charge, 0=0, 1=+, -1=-
    // io: order, 0=LO, 1=NLO
    
    void fdss_ (int &ih, int &ic, int &io, double &x, double& q2, double& u,
        double &ub, double &d, double &db, double &s, double &sb, double &c,
        double &b, double &g);

    extern struct{
        double fini;
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
        qs=1;
        //cerr << "Q=" << qs << " GeV out of Q^2 range [1,1e5] GeV^2. " << LINEINFO << endl;
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
		cerr << "WARNING! DSS HM does not give same result as ffgen???" << endl;
        return Evaluate(p, PIM, x, qs) + Evaluate(p, KM, x, qs);
    } else if (h==HP)
    {
		cerr << "WARNING! DSS HP does not give same result as ffgen???" << endl;
        return Evaluate(p, PIP, x, qs) + Evaluate(p, KP, x, qs);
    } 
    else if (h==PI)	// pi^+ + pi^-
    {    
		return Evaluate(p, PIP, x, qs) + Evaluate(p, PIM, x, qs);
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
    if (order==LO) return "DSS LO";
    return "DSS NLO";
}

DSS::DSS()
{
    initialized = false;
    fragini_.fini=0;
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
	
	SetOrder(NLO);
	double result, cor;
	result=0.1*Evaluate(U, PI0, 0.1, std::sqrt(10)); cor=0.64450;
	cout <<"z*D_{u->pi0}(Q^2=10GeV^2, x=0.1) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.2*Evaluate(D, PI, 0.2, std::sqrt(50)); cor=0.97394;
	cout <<"z*D_{d->pi+ pi-}(Q^2=50GeV^2, x=0.2) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.3*Evaluate(S, HM, 0.3, std::sqrt(20)); cor=0.64972;
	cout <<"z*D_{s->h-}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.3*Evaluate(G, PI0, 0.3, std::sqrt(20)); cor=0.37437;
	cout <<"z*D_{g->pi0}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
		
	cout <<"#LO:" << endl;
	SetOrder(LO);
	
		result=0.1*Evaluate(U, PI0, 0.1, std::sqrt(10)); cor=0.64300;
	cout <<"z*D_{u->pi0}(Q^2=10GeV^2, x=0.1) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.2*Evaluate(D, PI, 0.2, std::sqrt(50)); cor=1.0146;
	cout <<"z*D_{d->pi+ pi-}(Q^2=50GeV^2, x=0.2) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.3*Evaluate(S, HM, 0.3, std::sqrt(20)); cor=0.74499;
	cout <<"z*D_{s->h-}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	result=0.3*Evaluate(G, PI0, 0.3, std::sqrt(20)); cor=0.6663;
	cout <<"z*D_{g->pi0}(Q^2=20GeV^2, x=0.3) = " << result << " (correct " << cor << ")" <<  endl;
	if (std::abs(result-cor)/cor>0.01)
		cout << "TEST FAILED!!!" << endl;
	
	
	cout << "All tests done, if no errors were shown, all tests passed!" << endl;
	SetOrder(orig_order);
	
}
