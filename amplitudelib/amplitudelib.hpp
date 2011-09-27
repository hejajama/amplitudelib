/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "../tools/config.hpp"
#include "../tools/interpolation.hpp"
#include <vector>

class AmplitudeLib
{
    public:
        AmplitudeLib(std::string datafile);
        ~AmplitudeLib();

        // der w.r.t r der times.
        REAL N(REAL r, REAL y, int der=0);

        // S = 1-N
        REAL S(REAL r, REAL y, int der=0);
        
        // Amplitude in k-space, ft with 1/r^2 prefactor
        REAL N_k(REAL kt, REAL y);

        // Regular ft to k-space for S=1-N, normalization factor 1/(2\pi)^2
        REAL S_k(REAL kt, REAL y);

        // Virtual photon-proton cross sections, longitudinal and transverse
        // Notice that these are not normalized, as we don't integrate over
        // impact parameter
        // pol 0=longitudinal, 1=transverse
        REAL ProtonPhotonCrossSection(REAL Qsqr, REAL y, int pol);

        // Unintegrated gluon density
        REAL UGD(REAL k, REAL y, Interpolator* interp=NULL);

        // k_T factorization: d\sigma / (dyd^2p_T)
        // = const * 1/p_T^2 \int d^2 k_T/4 \alphas_(Q) \psi(|p_t+k_T|/2,x1)
        //   * \psi(|p_t-k_T|/2, x2)
        REAL dN_gluon_dyd2pt(REAL pt, REAL y, REAL sqrts);
        REAL dSigmady(REAL y, REAL sqrts);
        REAL dSigmady_mc(REAL y, REAL sqrts);

        // d ln N / d ln r^2
        REAL LogLogDerivative(REAL r, REAL y);

        // Saturation scale N(r, y) = Ns
        REAL SaturationScale(REAL y, REAL Ns);

        void InitializeInterpolation(REAL y);
        Interpolator* MakeInterpolator(REAL y);
    
        int YVals();
        int RPoints();
        REAL MinR();
        REAL MaxR();
        REAL MaxY();

        void SetOutOfRangeErrors(bool er);
        
        
    private:
        std::vector< std::vector<REAL> > n;
        std::vector<REAL> yvals;
        std::vector<REAL> lnrvals;
        std::vector<REAL> rvals;
        Interpolator *interpolator;

        REAL interpolator_y;
        REAL* tmprarray;
        REAL* tmpnarray;

        REAL minr;
        REAL rmultiplier;
        int rpoints;

        bool out_of_range_errors;  // don't print "out of range" errors
};

const int INTERPOLATION_POINTS = 8;
const REAL UGD_IR_CUTOFF=0.3;   // ugd(k<UGD_IR_CUTOFF)=0     BAD?????

#endif
