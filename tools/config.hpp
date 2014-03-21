/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
 */

#ifndef _CONFIG_HPP
#define _CONFIG_HPP 


#include <iostream>
using std::cout;
using std::cerr;
using std::cin;
using std::endl;
#include <string>
using std::string;
#include <cmath>

//typedef double REAL;
typedef unsigned int uint;


// Constants in Amplitude namespace, avoid overlapping with other programs
// using this class

namespace Amplitude
{

// Physical constants
//const REAL LAMBDAQCD2 = 0.21416*0.21416;   // GeV^2
const double LAMBDAQCD2 = 0.241*0.241;
const double LAMBDAQCD = 0.241;

const int Nf=3;
const int Nc=3;
const double Cf = (Nc*Nc-1.0)/(2.0*Nc);
const double ALPHA_e = 1.0/137.035999679; 
//const double e = sqrt(4.0*M_PI*ALPHA_e);

const double FMGEV = 5.068;



// Other constants

const double eps=0.000001;

}

// Inline functions

#define LINEINFO __FILE__ << ":" << __LINE__

inline double SQR(const double x) { return x*x; }

enum Parton
{
		    UVAL,DVAL,USEA,DSEA,U,D,S,C,B,G,    // Valence quarks, sea quarks,all quarks, gluons
            UBAR,DBAR,SBAR, LIGHT			// LIGHT = all light quarks
};

enum Order
{
    LO,
    NLO
};

#endif
