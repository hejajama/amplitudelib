#ifndef _SINGLE_INCLUSIVE_HPP
#define _SINGLE_INCLUSIVE_HPP

/*
 * Single inclusive cross section
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2015
 */
 

#include "amplitudelib.hpp"
#include "../tools/config.hpp"
#include <vector>

/**
 * Calculate single inclusive yields
 *
 * Calculate single particle production using the hybrid formalism
 * or the kt factorization.
 *
 */
class SingleInclusive
{
    public:

        SingleInclusive(AmplitudeLib* N_);
        
        /**
         * Single inclusive hadron yield
         *
         * Compute single inclusive hadron yield dN/(d^2 pt dy) using the
         * hybrid formalism
         * @param y Rapidity of the produced particle
         * @param pt Transverse momentum of the produced particle
         * @param sqrts Center-of-mass energy
         * @param fragfun Pointer to the fragmentation function class
         * @param pdf Pointer to the PDF class
         * @param deuteron If true, probe is deuteron=proton+neutron
         * @param final Produced particle (default: pi0)
         * @param scale Scale at which PDF and FF are evaluated, default: transverse momentum
         * @see PDF
         * @see FragmentationFunction
         */
        double dHadronMultiplicity_dyd2pt(double y, double pt, double sqrts,
            FragmentationFunction *fragfun, PDF* pdf, Amplitude::Hadron final=Amplitude::PI0,
            bool deuteron=false,  double scale=-1 );

        

        /**
         * Single inclusive hadron yield, kt-factorization
         *
         * Compute single inclusive hadron yield dN/(d^2 pt dy) using the
         * kt-factorization. Ref e.g. hep-ph/0111362 (40)
         * @param y Rapidity of the produced particle
         * @param pt Transverse momentum of the produced particle
         * @param sqrts Center-of-mass energy
         * @param fragfun Pointer to the fragmentation function class
         * @param final Produced particle (default: pi0)
         * @param N2 Dipole amplitude for the probe (the saved dipole amplitude is for the target), default: use the same for both
         * @see PDF
         * @see FragmentationFunction
         */
		double dHadronMultiplicity_dyd2pt_ktfact(double y, double pt, double sqrts, FragmentationFunction* fragfun, Amplitude::Hadron final, AmplitudeLib* N2=NULL );

        /**
         * Parton level single inclusive yield using kt-factorization
         *
         * Compute single parton yield using kt-factorization
         * @param y Rapidity of the produced particle
         * @param pt Transverse momentum of the produced particle
         * @param sqrts Center-of-mass energy
         * @param N2 Dipole amplitude for the probe (the saved dipole amplitude is for the target), default: use the same for both
         * @param scale Momentum scale [GeV^2] at which the strong coupling constant is evaluated, default: parton pt^2
         */
		double dHadronMultiplicity_dyd2pt_ktfact_parton(double y, double pt, double sqrts, AmplitudeLib* N2=NULL, double scale=-1 );

        /**
         * Set partons included in hybrid formalism calclulation
         *
         * @param p: List of partons. Possible values: LIGHT, C, B, G. Default: LIGHT and G
         */
        void SetPartons(std::vector<Amplitude::Parton> p);

        /**
         * Get list of partons taken into account in hybrid formalism
         */
        std::vector<Amplitude::Parton> &Partons();

        /*
         * TODO:
        
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

		*/
    private:
        AmplitudeLib* N;
        QCD qcd;
        std::vector<Amplitude::Parton> partons;            //! Partons included in hybrid formalism

};

#endif
		
