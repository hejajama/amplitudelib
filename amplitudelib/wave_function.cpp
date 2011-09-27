/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010
 */
 
 
#include "wave_function.hpp"
#include <iostream>

WaveFunction::WaveFunction()
{

}


double WaveFunction::PsiSqr_tot(double Qsqr, double r, double z)
{
    return PsiSqr_T(Qsqr,r,z)+PsiSqr_L(Qsqr,r,z);
}

double WaveFunction::PsiSqr_tot_intz(double Qsqr, double r)
{
    return PsiSqr_T_intz(Qsqr,r)+PsiSqr_L_intz(Qsqr,r);
}

