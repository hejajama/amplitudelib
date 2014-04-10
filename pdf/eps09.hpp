#ifndef _EPS09_H
#define _EPS09_H

/*
 * Wrapper calss for EPS09 nuclear parton distribution functions 
 * Real code is in file EPS09.f downloaded from
 * https://www.jyu.fi/fysiikka/en/research/highenergy/urhic/eps09
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013
 */

#include "../tools/config.hpp"
#include "pdf.hpp"
#include "cteq.hpp"
#include <string>

// Nuclear mass number A is set using PDF->SetA() saved in 
// protected variable int A
class EPS09 : public PDF
{
    public:
        double xq(double x, double q, Amplitude::Parton p);    // return x*q(x,q)
        void Initialize(int param=-1);
        std::string GetString();
        
        void SetOrder(Amplitude::Order o);

        double MinQ();
        double MaxQ();
        
        EPS09();
        int SetA(int A_);
        
        void Test();
    private:
		CTEQ cteq;	// eps09 only gives nuclear modification factor, it is multiplied by x*f(x) from cteq


};

// The following functions are implemented in CT10Pdf.f file
extern "C" 
{
        void eps09_(int& order, int& set, int& A, double& x, double &Q, double& ruv,
            double& rdv, double& ru, double& rd, double& rs, double& rc, double& rb, double& rg);

}

#endif
