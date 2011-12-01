#ifndef _HKNS_H
#define _HKNS_H

/*
 * HKNS fragmentation function
 * Uses fragmentation_hknsff07.f downloaded from
 * http://research.kek.jp/people/kumanos/ffs.html
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

class HKNS : public FragmentationFunction
{
    public:
        REAL Evaluate(Parton p, Hadron h, REAL x, REAL qs);
        HKNS();
        std::string GetString();
    private:


};


#endif
