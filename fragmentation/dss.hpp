#ifndef _DSS_H
#define _DSS_H

/*
 * DSS fragmentation function
 * Uses fragmentation_dss.f, ref hep-ph/0703242
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <string>
#include "../tools/config.hpp"
#include "fragmentation.hpp"

class DSS : public FragmentationFunction
{
    public:
        double Evaluate(Parton p, Hadron h, double x, double qs);
        DSS();
        std::string GetString();
        void SetOrder(Order order_);
        
        void Test();
    private:
        bool initialized;

};


#endif
