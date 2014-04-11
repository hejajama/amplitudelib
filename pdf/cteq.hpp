#ifndef _CTEQ_H
#define _CTEQ_H

/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include "../tools/config.hpp"
#include "pdf.hpp"
#include <string>

/**
 * Wrapper calss for CTEQ-TEA parton distribution functions CT12
 * Real code is in file CT12Pdf.f downloaded from
 *
 * For LO distributions Cteq6pdf.f is used.
 * 
 * http://hep.pa.msu.edu/cteq/public/index.html
 *
 */
class CTEQ : public PDF
{
    public:
        double xq(double x, double q, Amplitude::Parton p);    // return x*q(x,q)

        std::string GetString();

        /**
         * Select between NLO and LO set
         */
        void SetOrder(Amplitude::Order o);

        double MinQ();
        double MaxQ();
        
        CTEQ();
        
        void Test();


};

extern "C"
{
    
    void setctq6_(int& iset_);
    double ctq6pdf_(int& iparton, double& x, double& q);
    
    void setct12_(char* file);		// file is 40-char string!
    double ct12pdf_(int& iparton, double& x, double& q);
}

#endif
