#ifndef _PDF_H
#define _PDF_H

/*
 * Virtual class to hide different parton distribution functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>
 */

#include "../tools/config.hpp"
#include <string>

class PDF
{
    public:
		PDF();
        ~PDF();
        virtual double xq(double x, double q, Parton p)=0;    // return x*q(x,q)
        // x_1*x_2*f(x_1,x_2) dpfd with kinematical constraint x_1+x_2 < 1
        double Dpdf(double x1, double x2, double q, Parton p1, Parton p2);	
        virtual void Initialize(int param=-1);
        virtual std::string GetString();
        virtual double MinX();
        virtual double MinQ();   // GeV
        virtual double MaxX();
        virtual double MaxQ();

        void PlotPdf(double Q);
	protected:
		bool initialized;
};

#endif
