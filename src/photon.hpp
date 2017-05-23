#ifndef _SINGLE_INCLUSIVE_HPP
#define _SINGLE_INCLUSIVE_HPP

/*
 * Inclusive isolated photon calculation
 * 
 * Heikki Mäntysaari <mantysaari@bnl.gov>, 2017
 */

// General kinematics:
// Photon momenta k, rapidity y_k
// Quark momenta l, rapidity y_l
// phi: angle between k and l


#include "../amplitudelib/amplitudelib.hpp"
#include "../tools/config.hpp"
#include <vector>

using namespace Amplitude;
/**
 * Calculate inclusive isolated photon
 *
 */
class IsolatedPhoton
{
public:
    
    IsolatedPhoton(AmplitudeLib* N_) { N=N_; sqrts=200; partons.push_back(U); partons.push_back(D); partons.push_back(S); partons.push_back(C); };

    void SetIsolationCut(double R) { isolation_cut = R; };
    double DifferentialPhotonCrossSection(double k, double y_k, double l, double y_l, double phi);
    double PhotonCrossSection(double k, double y_k); // d sigma / d^2 k dy_k
    
    AmplitudeLib* GetDipole() { return N; }
    std::vector<Amplitude::Parton> GetPartons() { return partons; }
    
private:
    AmplitudeLib* N;
        QCD qcd;
        double sqrts;
        std::vector<Amplitude::Parton> partons;            //! Partons included
    
        double isolation_cut;
    
};

#endif

