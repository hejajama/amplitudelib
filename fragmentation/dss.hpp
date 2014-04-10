#ifndef _DSS_H
#define _DSS_H

/*
 * DSS fragmentation function
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

/**
 * DSS fragmentation function
 *
 * Ref. hep-ph/0703242
 */

class DSS : public FragmentationFunction
{
    public:
        /**
         * Evaluate FF
         *
         * Evaluate DSS fragmentation function D(z)
         * @param p Parton
         * @param h Hadron
         * @param x Longitudinal momentum fraction
         * @param qs Scale [GeV]
         */
        double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double x, double qs);
        DSS();
        std::string GetString();
        void SetOrder(Amplitude::Order order_);
        
        void Test();
    private:
        bool initialized;

};


#endif
