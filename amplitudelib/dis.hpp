#ifndef _DIS_HPP
#define _DIS_HPP

#include "amplitudelib.hpp"

/**
 * Deep inelastic scattering
 *
 * Computes proton structure functions (F2, FL) and the reduced
 * cross section.
 */
class DIS
{
    public:
        /**
         * Constructor, initialize calculation
         *
         * Initialize DIS calculations. Requires a pointer to the AmplitudeLib object
         */
         DIS(AmplitudeLib* amp);

        // Virtual photon-proton cross sections, longitudinal and transverse
        // Notice that these are not normalized, as we don't integrate over
        // impact parameter
        /**
         * Virtual photon-proton cross section
         *
         * Computes virtual photon-proton total cross section at given kinematics and photon polarization
         * @param qsqr Photon virtuality (GeV^2)
         * @param xbj Bjorken-x
         * @param pol Photon polarization
         * @param Parton quark flavor(s) taken into account (default: U,D,S)
         * @param mass Quark mass, if not specified, default value from VirtualPhoton class is used
         * @see VirtualPhoton
         */ 
        double ProtonPhotonCrossSection(double qsqr, double xbj, Amplitude::Polarization pol, Amplitude::Parton=LIGHT, double mass=-1);

        /**
         * Structure function F2
         *
         * Compute structure function F2
         * @param qsqr Photon virtuality (GeV^2)
         * @param xbj Bjorken-x
         * @param Parton quark flavor(s) taken into account (default: U,D,S)
         * @param mass Quark mass, if not specified, default value from VirtualPhoton class is used
         */
        double F2(double qsqr, double xbj, Parton=LIGHT, double mass=-1);

        /**
         * Structure function FL
         *
         * Compute structure function FL. Same parameters as in F2()
         * @see F2
         */
        double FL(double qsqr, double xbj, Parton=LIGHT, double mass=-1);

        /**
         * Compute reduced cross section
         *
         * Computes the reduced photon-proton cross section
         * sigma_r = F2 - y^2 / ( 1.0 + (1.0-y)^2 ) * FL;
         * @param qsqr Photon virtuality (GeV^2)
         * @param xbj Bjorken-x
         * @param sqrts Center-of-mass energy sqrt(s)
         * @param Parton quark flavor(s) taken into account (default: U,D,S)
         * @param mass Quark mass, if not specified, default value from VirtualPhoton class is used
         * @see F2
         * @see FL
         */
        double ReducedCrossSection(double qsqr, double xbj, double sqrts, Parton=LIGHT, double mass=-1);

    private:
        AmplitudeLib* N;
};

#endif
