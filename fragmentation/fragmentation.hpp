#ifndef _FRAGMENTATION_H
#define _FRAGMENTATION_H

/*
 * Virtual class to hide different fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "../tools/config.hpp"


class FragmentationFunction
{
    public:
        // D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
        virtual double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double x, double qs);
        FragmentationFunction();
        virtual std::string GetString();
        Amplitude::Order GetOrder();
        virtual void SetOrder(Amplitude::Order order_);
        
        virtual void Test();
    protected:
        Amplitude::Order order;  // 0: LO, 1: NLO

};


std::string ParticleStr(Amplitude::Hadron h);

#endif
