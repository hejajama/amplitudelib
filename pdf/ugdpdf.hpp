#ifndef _UGDPDF_H
#define _UGDPDF_H

/*
 * PDF from UGD
 * Computes UGD from dipole amplitude and integrates it over k_T
 * to get PDF
 * 
 * Note: requires pointer to AmplitudeLib* and the scaling
 * factor \sigma_0/2
 * 
 * NOTE: Deletes the AmplitudeLib object in destructor! So assumes
 * that AmplitudeLib* N is not used outside this class!
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#include "../tools/config.hpp"
#include "pdf.hpp"
#include "../amplitudelib/amplitudelib.hpp"
#include <string>

class UGDPDF : public PDF
{
    public:
        double xq(double x, double q, Parton p);    // return x*q(x,q)
        //void Initialize(int param=-1);
        std::string GetString();
        
        double MaxX();
        double MinX();
        double MaxQ();
        double MinQ();
        
        UGDPDF();
        UGDPDF(AmplitudeLib* N, double sigma02_=1.0);
        ~UGDPDF();
        
        //void Test();
    private:
		AmplitudeLib* N;
		double sigma02 ;


};


#endif

