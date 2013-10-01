/*
 * AmplitudeLib for the dipole amplitude
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
 */

#ifndef _AMPLITUDELIB_H
#define _AMPLITUDELIB_H

#include "../tools/config.hpp"
#include "../tools/interpolation.hpp"
#include "../fragmentation/fragmentation.hpp"
#include "../pdf/pdf.hpp"
#include <vector>
#include <string>

enum Polarization
{
	L,
	T
};

// Running/fixed coupling
enum RUNNING_ALPHAS
{
	RUNNING,
	FIXED
};

enum FT_METHOD
{
	GSL,
	ACC_SERIES   // fourier/fourier.c
};

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

        // Regular ft to k-space for S^p=(1-N)^p, NO normalization factor
        // \int e^(ik.r) S(r)
        // if adjoint=true, then N is in adjoint representation
        double S_k(double kt, double y, bool adjoint=false, double pow=1.0);
        
        

        // Amplitude in adjoint representation
        double N_A(double r, double y, int der=0);
        
        void SetX0(double x0_);	// override x0


        // Virtual photon-proton cross sections, longitudinal and transverse
        // Notice that these are not normalized, as we don't integrate over
        // impact parameter
        
        double ProtonPhotonCrossSection(double Qsqr, double y, Polarization pol, Parton=LIGHT); // default: sum over u,d,s quarks
        double F2(double Qsqr, double y, Parton=LIGHT);
        double FL(double Qsqr, double y, Parton=LIGHT);
        double ReducedCrossSection(double qsqr, double y, double sqrts, Parton=LIGHT);


		//////////////// SINGLE INCLUSIVE, hybrid formalism
		//////// File xs.cpp

        // Differential forward hadron production multiplicity
        // dN_h / (dy_h d^2 p_T)
        double dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0, double scale=-1 );
        // Parton level
        double dHadronMultiplicity_dyd2pt_parton(double y, double pt, double sqrts,
			PDF* pdf, bool deuteron, double scale=-1);
        // Integrated over rapidity and pt range
        double HadronMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0 );
        // Average hadron multiplicity in rapidity range
        double AverageHadronMultiplicity(double miny, double maxy, double pt, double sqrts, 
            FragmentationFunction *fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0 );
        // Douple parton scattering at fixed pt, y
        double DPS(double y1, double y2, double pt1, double pt2, double sqrts,
              FragmentationFunction* fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0, char dps_mode='c');
        // Parton level DPS, includes PDF and sum over quarks/gluons, doesn't include fragmentation
        double DPS_partonlevel(double y1, double y2, double pt1, double pt2, double sqrts,
              PDF* pdf, bool deuteron=false, char dps_mode='c', double scale=-1);
        double DPSMultiplicity(double miny, double maxy, double minpt, double maxpt, double sqrts,
			FragmentationFunction* fragfun, PDF* pdf, bool deuteron=false, Hadron final=PI0, char dps_mode='c');

		
		/////////////// UGD  (file ugd.cpp)
		
		// KMR UGD C_F/(8\pi^3) S_T/\alpha_s(q) q^4 S_k(q)
		// returns without factor S_T unless it is specified
		double UGD(double q, double y, double scale_=-1, double S_T=1.0);
		// Integrated GD from UGD
		double xg(double x, double q);
		
		// k_T factorization dN/d^2pt dy, ref e.g. hep-ph/0111362 (40)
		double dHadronMultiplicity_dyd2pt_ktfact(double y, double pt, double sqrts, FragmentationFunction* fragfun, Hadron final, AmplitudeLib* N2=NULL );
		double dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2=NULL );
		
		
		///////


        // d ln N / d ln r^2
        double LogLogDerivative(double r, double y);

        // Saturation scale N(r, y) = Ns
        // Returns r
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
        
        double Sigma02();
        void SetSigma02(double s_);

        double X0();

        bool SetOutOfRangeErrors(bool er);
        
        std::string GetString();
        
        
        double Alphas(double qsqr);
        
        void SetRunningCoupling(RUNNING_ALPHAS as_);
        RUNNING_ALPHAS GetRunningCoupling();
        
        std::string Version();
        
        void SetFTMethod(FT_METHOD f);
        FT_METHOD GetFTMethod();		

        
        
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
        
        double sigma02;		// sigma_0 / 2,  ktfactorization hadronprod results are multiplied by this

        bool out_of_range_errors;  // don't print "out of range" errors
        
        std::string info_string;
        
        RUNNING_ALPHAS as;
        FT_METHOD ft;  // ACC SERIES: use j0_transfer from fourier/fourier.c,	
					// it should be faster but sometimes it is much slower!!
        
        
};

const int INTERPOLATION_POINTS = 12;
const double UGD_IR_CUTOFF=0.3;   // ugd(k<UGD_IR_CUTOFF)=0     BAD?????

const int FOURIER_ZEROS=2000;   // How many zeros of the Bessel functions is
                    // used when Fourier transforming



const std::string AMPLITUDELIB_VERSION = "1.1-dev 2013-06-xx";
#endif
