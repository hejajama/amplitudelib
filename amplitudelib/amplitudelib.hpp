/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "../tools/config.hpp"
#include "../tools/interpolation.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../pdf/pdf.hpp"
#include <vector>

class AmplitudeLib
{
    public:
        AmplitudeLib(std::string datafile, bool kspace_=false);
        ~AmplitudeLib();

        // der w.r.t r der times.
        double N(double r, double y, int der=0, bool bspline=false);

        // S = 1-N
        double S(double r, double y, int der=0);
        
        // Amplitude in k-space, ft with 1/r^2 prefactor
        double N_k(double kt, double y);

        // Amplitude from k space to x space
        double N_k_to_x(double x, double y);

        // Regular ft to k-space for S=1-N, NO normalization factor
        // \int e^(ik.r) S(r)
        // if adjoint=true, then N is in adjoint representation
        double S_k(double kt, double y, bool adjoint=false);

        // Amplitude in adjoint representation
        double N_A(double r, double y, int der=0);


        // Virtual photon-proton cross sections, longitudinal and transverse
        // Notice that these are not normalized, as we don't integrate over
        // impact parameter
        // pol 0=longitudinal, 1=transverse
        double ProtonPhotonCrossSection(double Qsqr, double y, int pol);

        // Differential forward hadron production multiplicity
        // dN_h / (dy_h d^2 p_T)
        double dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0 );
        // Integrated over rapidity and pt range
        double HadronMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0 );
        // Douple parton scattering, integrate products of two dHadronMultiplicity_dyd2pt
        // over kinematical regions
        // pt2<pt1
        double DPS(double miny, double maxy, double minpt1, double minpt2, double sqrts,
            FragmentationFunction* fragfun, bool deuteron=false, Hadron final=PI0);

        // Unintegrated gluon density
        double UGD(double k, double y, Interpolator* interp=NULL);

        // k_T factorization: d\sigma / (dyd^2p_T)
        // = const * 1/p_T^2 \int d^2 k_T/4 \alphas_(Q) \psi(|p_t+k_T|/2,x1)
        //   * \psi(|p_t-k_T|/2, x2)
        double dN_gluon_dyd2pt(double pt, double y, double sqrts);
        double dSigmady(double y, double sqrts);
        double dSigmady_mc(double y, double sqrts);

        // d ln N / d ln r^2
        double LogLogDerivative(double r, double y);

        // Saturation scale N(r, y) = Ns
        double SaturationScale(double y, double Ns);

        void InitializeInterpolation(double y, bool bspline=false);
        bool InterpolatorInitialized(double y); // Returns whether or not the interpolator
                                    // is initialized
        Interpolator* MakeInterpolator(double y);
    
        int YVals();
        int RPoints();
        double MinR();
        double MaxR();
        double MaxY();

        double X0();

        bool SetOutOfRangeErrors(bool er);
        
        
    private:
        // [yind][r/kind]
        std::vector< std::vector<double> > n;
        std::vector<double> yvals;
        std::vector<double> lnrvals;
        std::vector<double> rvals;
        Interpolator *interpolator;

        bool kspace;    // true if data is in kspace

        double interpolator_y;
        double* tmprarray;
        double* tmpnarray;

        double minr;
        double rmultiplier;
        
        double maxr_interpolate;  // When interpolating in coordinate space,
                    // force N(r>maxr_interpolate)=1 (avoid osciallatory 
                    // artifacts from interpolation code)
        int rpoints;

        double x0;

        bool out_of_range_errors;  // don't print "out of range" errors
};

const int INTERPOLATION_POINTS = 12;
const double UGD_IR_CUTOFF=0.3;   // ugd(k<UGD_IR_CUTOFF)=0     BAD?????

const int FOURIER_ZEROS=700;   // How many zeros of the Bessel functions is
                    // used when Fourier transforming
#endif
