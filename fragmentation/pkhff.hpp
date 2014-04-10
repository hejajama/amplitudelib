#ifndef _PKHFF_H
#define _PKHFF_H

/*
 * PKHFF fragmentation function
 * Supports only charged hadrons
 * Uses pkhff.f downloaded from
 * http://www2.pv.infn.it/~radici/FFdatabase/
 * Ref Phys. Rev. D 62 (2000) 054001
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

class PKHFF : public FragmentationFunction
{
    public:
        double Evaluate(Amplitude::Parton p, Amplitude::Hadron h, double x, double qs);
        PKHFF();
        std::string GetString();
    private:
        int ifini;  // whether or not PKHFF is initialized, should be 0 initially,
                // fortran code updates

};


#endif
