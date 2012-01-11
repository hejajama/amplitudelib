#ifndef _FRAGMENTATION_H
#define _FRAGMENTATION_H

/*
 * Virtual class to hide different fragmentation functions
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../tools/config.hpp"

enum Hadron
{
    PI,   // pi+, pi-
    PIP,    // pi^+
    PIM,    // pi^-
    K,      // k+, k-
    KP,     // k^+
    KM,     // k^-
    K0,      // K^0, \bar K^0
    P,      // P, \bar P
    PP,     // p^-
    PM,     // p^+
    PI0,    // \pi^0
    NE,      // N, \bar N
    H,       // h^+ + h^- sum of charged hadrons
    HP,     // h^+
    HM      // h^-
};

enum Order
{
    LO,
    NLO
};

class FragmentationFunction
{
    public:
        // D_{p->h}, x: long. mom. fraction, qs: scale (GeV)
        virtual REAL Evaluate(Parton p, Hadron h, REAL x, REAL qs);
        FragmentationFunction();
        virtual std::string GetString();
        Order GetOrder();
        void SetOrder(Order order_);
    protected:
        Order order;  // 0: LO, 1: NLO

};


#endif