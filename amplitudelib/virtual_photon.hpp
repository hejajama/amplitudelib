#ifndef VirtualPhoton_H
#define VirtualPhoton_H

/*
 * Virtual photon-q\bar q overlap
 * Ref: Kowalski, Motyka and Watt, see arXiv: hep-ph/0606272v2
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2011
 */


#include "wave_function.hpp"
#include <iostream>
#include <string>

// f = quark flavor, 0=u, 1=d, 2=s

class VirtualPhoton : public WaveFunction {
    public:
        VirtualPhoton();
        
        // Overlap wave functions
        double PsiSqr_T(double Qsqr, double r, double z);
        double PsiSqr_L(double Qsqr, double r, double z);
        
        // Overlap wave functions integrated over z=[0,1]
        double PsiSqr_T_intz(double Qsqr, double r);
        double PsiSqr_L_intz(double Qsqr, double r);
        
        
        std::string GetParamString();
        
        
    private:
        // Parameters
        double e_f[3];   // Quark charges
        double m_f[3];   // Quark masses (GeV)

        double Epsilon(double Qsqr, double z, int f);
        
        
};

std::ostream& operator<<(std::ostream& os, VirtualPhoton& ic);

#endif
