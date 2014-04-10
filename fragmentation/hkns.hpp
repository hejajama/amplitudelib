#ifndef _HKNS_H
#define _HKNS_H

/*
 * HKNS fragmentation function
 * Uses fragmentation_hknsff07.f downloaded from
 * http://research.kek.jp/people/kumanos/ffs.html
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

class HKNS : public FragmentationFunction
{
    public:
        double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double x, double qs);
        HKNS();
        std::string GetString();
    private:


};


#endif
