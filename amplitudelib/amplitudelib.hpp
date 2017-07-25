/*
 * AmplitudeLib2 for the dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2014
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "../tools/config.hpp"
#include "../tools/interpolation.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../pdf/pdf.hpp"
#include "qcd.hpp"
#include <vector>
#include <string>





/**
 * Amplitude class
 *
 * Stores the loaded dipole amplitude (loaded by DataFile class), and
 * evaluates the amplitude at given r,xbj.
 * Also computes Fourier transfers between transverse coordinate and
 * momentum spaces.
 */
class AmplitudeLib
{
    public:
        /**
         * Loader.
         * 
         * Constructor, loads the solution for the BK equation from
         * the given file and stores the result in memory.
         * If the second optional argument is true, then it is
         * assumed that the dipole amplitude is solved in momentum
         * space
         */
        AmplitudeLib(std::string datafile, bool kspace_=false);
    
        /* Load amplitude from the give data tables
         *
         * data[i][j] is dipole amplitude at rapidity yvals[i] for dipole size rvals[i]
         */
        AmplitudeLib(std::vector< std::vector< double > > &data, std::vector<double> &yvals_, std::vector<double> &rvals_);
    
        ~AmplitudeLib();

        /**
         * Dipole amplitude
         *
         * Evaluates the dipole amplitude at given dipole size r and at
         * given Bjorken x (xbj). In momentum space the first argument is
         * transverse momentum.
         */
        double N(double r, double xbj);

        /**
         * Dipole amplitude in adjoint representation
         *
         * Evaluates the dipole amplitude in adjoint representation,
         * 2N(r)-N(r)^2, at given dipole size r and at
         * given Bjorken x (xbj). In momentum space the first argument is
         * transverse momentum.
         */
        double N_A(double r, double xbj);

        /**
         * Scattering matrix
         *
         * Evaluate the scattering matrix 1-N (forced between 0 and 1)
         * at given dipole size and Bjorken-x.
         */
         double S(double r, double xbj);
        
        /**
         * Weizsäcker-Williams (WW) gluon distribution
         *
         * Computes the Weizsäcker-Williams (WW) gluon distribution
         * (in momentum space), defiend as
         *   WW_ugd = \int d^2 r/(2 pi r^2) exp(ik.r)N(r)
         *       = \int dr/r BesselJ[0,k*r] * N(r)
         * @param kt transverse momentum
         * @param xbj Bjorken-x
         */
        double WW_UGD(double kt, double xbj);

        /**
         * Fourier transfer to coordinate space
         *
         * Computes the dipole amplitude in coordinate space, defiend as
         * 
         *  N(r) = r^2 \int d^2 k/(2 pi) exp(-ik.r) N(k)
         *       = r^2 \int dk k BesselJ[0,k*r] * N(k)
         *
         * Note: requires that the loaded bk solution is in momentum space
         * @param r dipole size
         * @param xbj Bjorken-x
         */
        double N_k_to_x(double r, double xbj);

        /**
         * Scattering matrix in momentum space
         *
         * Compute (by Fourier transform) the scattering matrix in
         * momentum space
         * 
         *  S_k = \int e^(ik.r) S(r)^pow
         *
         * If representation is set to ADJOINT, computes using the
         * adjoint representation dipole amplitude.
         *
         * By default p=1, but user can set different pow, in which
         * case fourier transfer of (1-N(r))^pow is computed. 
         */
        double S_k(double kt, double xbj, Amplitude::Representation=Amplitude::FUNDAMENTAL, double pow=1.0);
        

        /**
         * Dipole unintegrated gluon distribution
         *
         * Defined as
         * 
		 * C_F/(8 pi^3) S_T/alpha_s(q) q^4 S_k(q)
         *
         * S_k is computed in adjoint representation. 
		 * Computes without factor S_T unless it is specified
         * @param q Scale (GeV)
         * @param as_scale scale at which the running coupling is evaluated, default: q
         * @param S_T Transverse size of the target = normalization factor
         * @see S_k
         * 
         */
		double Dipole_UGD(double q, double xbj, double as_scale_=-1, double S_T=1.0);

        /**
         * Gluon distribution from UGD
         *
         * Integrated gluon distribution computed from the unintegrated gluon distribution
         * Computed as
         *
         * int_0^Q dq^2/q^2 UGD(q)
         * @param x Bjorken-x
         * @param q Scale
         */
		double xg(double x, double q);

        /**
         * Initialize interpolation for certain Bjorken-x
         */
        void InitializeInterpolation(double xbj);

        /**
         * Test if interpolator is initialized at given Bjorken-x
         */
        bool InterpolatorInitialized(double xbj);

        /**
         * Create interpolator at given Bjorken-x and return it
         */
        Interpolator* MakeInterpolator(double xbj);

        // d ln N / d ln r^2
        double LogLogDerivative(double r, double y);

        /**
         * Saturation scale
         *
         * Solve saturation scale defined as
         * 
         * N(r, xbj)=Ns
         * 
         * @return Dipole size r that satisfies the satscale condition (momentum space: TODO)
         */
        double SaturationScale(double xbj, double Ns);

        
        /**
         * Number or rapidity values in the BK solution
         */
        int YPoints();

        /**
         * Return ith rapidity value
         *
         * Can be used to avoid interpolation uncertainties 
         */
        double YValue(int yind);
         

        /**
         * Number of dipole sizes in the BK solution
         */
        int RPoints();

        /**
         * Minimum dipole size
         */
        double MinR();

        /**
         * Maximum dipole size
         */
        double MaxR();
        double MaxY();
        
        /**
         * Initial Bjorken-x in the BK solution
         *
         * Can we owerwritten using SetX0() method
         * @see SetX0
         */
        double X0();

        /**
         * Set initial Bjorken-x
         *
         * Override the Bjoken-x at initial condition for the BK solution
         * (x0).
         */
        void SetX0(double x0_);

        /**
         * Specify whether an error message to stderr is printed for too
         * small/alrge dipoles or small/large xbj values
         *
         * @return previous setting
         */

        bool SetOutOfRangeErrors(bool er);

        /**
         * Returns an info string describing the BK solution and AmplitudeLib version
         */
        std::string GetString();

        /**
         * Return the Fourier transfer method used
         */
        Amplitude::FT_Method GetFTMethod();

        /**
         * Set Fourier transfer method
         * @see FT_Method
         */
        void SetFTMethod(Amplitude::FT_Method f);

        /**
         * Return version string
         */
        std::string Version();

        QCD& GetQCD();

    private:
        // [yind][r/kind]
        std::vector< std::vector<double> > n;
        std::vector<double> yvals;
        std::vector<double> lnrvals;
        std::vector<double> rvals;
        Interpolator *interpolator;

        bool kspace;    //! true if data is in kspace

        double interpolator_xbj;  //! xbj at which the interpolator is initialized
        double* tmprarray;
        double* tmpnarray;

        double minr;
        double rmultiplier;

        /**
         * Max dipole size for interpolation
         *
         * When interpolating in coordinate space,
         * force N(r>maxr_interpolate)=1 (avoid osciallatory 
         * artifacts from interpolation code)
         */
        double maxr_interpolate;
        
        int rpoints;    //! Number of dipole sizes

        double x0;      //! Initial Bjorken-x
        
        

        bool out_of_range_errors;  //! If true, don't print "out of range" errors
        
        std::string info_string;
 
        std::string datafilename;   //! Name of the file where the amplitude is read, for error messages

        QCD qcd;
        
        Amplitude::FT_Method ft;  // ACC SERIES: use j0_transfer from fourier/fourier.c,	
					// it should be faster but sometimes it is much slower!!
        
        
};




const int INTERPOLATION_POINTS = 12;

const int FOURIER_ZEROS=1000;   // How many zeros of the Bessel functions is
                    // used when Fourier transforming



const std::string AMPLITUDELIB_VERSION = "2.1 2017-07-25";

const Amplitude::FT_Method DEFAULT_FT_METHOD = Amplitude::ACC_SERIES;
#endif
