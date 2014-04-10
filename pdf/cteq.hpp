#ifndef _CTEQ_H
#define _CTEQ_H

/*
 * Wrapper calss for CTEQ-TEA parton distribution functions CT10
 * Real code is in file CT10Pdf.f downloaded from
 * http://hep.pa.msu.edu/cteq/public/index.html
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include "../tools/config.hpp"
#include "pdf.hpp"
#include <string>

class CTEQ : public PDF
{
    public:
        double xq(double x, double q, Amplitude::Parton p);    // return x*q(x,q)
        void Initialize(int param=-1);
        std::string GetString();
        
        void SetOrder(Amplitude::Order o);

        double MinQ();
        double MaxQ();
        
        CTEQ();
        
        void Test();


};

// The following functions are implemented in CT10Pdf.f file
extern "C"
{
    //void setct10_(int& iset_);
    //double ct10pdf_(int& iparton, double& x, double& q);
    
    void setctq6_(int& iset_);
    double ctq6pdf_(int& iparton, double& x, double& q);
    
    void setct12_(char* file);		// file is 40-char string!
    double ct12pdf_(int& iparton, double& x, double& q);
}

#endif
