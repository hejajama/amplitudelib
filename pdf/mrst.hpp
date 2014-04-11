#ifndef _MRST_H
#define _MRST_H

/*
 * AmplitudeLib
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012-2014
 */


#include "pdf.hpp"
#include "../tools/config.hpp"
#include "mrst99.h"
#include <string>

/**
 * Wrapper class to use mrst99 NLO PDF,
 * actual code is in files mrst99.{cpp,h} downloaded from
 * http://durpdg.dur.ac.uk/HEPDATA/PDF
 */
class MRST : public PDF
{
    public:
        ~MRST();
        double xq(double x, double q, Amplitude::Parton p);    // return x*q(x,q)

        /**
         * Initialize, param is given to the MRST code, and it
         * specifies the pdf set used. If param is not specified,
         * the default set is used. Sets are defined in file
         * mrst99.h
         */
        void Initialize(int param=-1);
        std::string GetString();
    private:
        c_mrst* mrst;
        int set;
};

#endif
