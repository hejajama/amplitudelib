#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H

/*
 * General class for wave functions
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2010-2011
 */

#include <string>

class WaveFunction{
    public:
        WaveFunction();
        virtual double PsiSqr_T(double Qsqr, double r, double z) = 0;
        virtual double PsiSqr_L(double Qsqr, double r, double z) = 0;
        virtual double PsiSqr_T_intz(double Qsqr, double r) = 0;
        virtual double PsiSqr_L_intz(double Qsqr, double r) = 0;
        virtual std::string GetParamString()=0;
        double PsiSqr_tot(double Qsqr, double r, double z);
        double PsiSqr_tot_intz(double Qsqr, double r);
    private:
        int mode;   // What to return when PsiSqr_intz
};

#endif  // WAVE_FUNCTION_H
