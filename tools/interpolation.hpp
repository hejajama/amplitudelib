/*
 * General purpose interpolation class
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2012
 */

#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

/*
 * Interpolates given data using spline (goes trough every data point)
 * or bspline (=noisy data)
 * Uses GSL
 *
 * Xdata and ydata pointers are saved, but not used for interpolation purposes
 * If user frees the allocated memory, one should be sure that these pointers
 * are not asked from this class!
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include "config.hpp"
#include "tools.hpp"

enum INTERPOLATION_METHOD {
    INTERPOLATE_SPLINE,
    INTERPOLATE_BSPLINE
};

class Interpolator
{
    public:
        Interpolator(REAL* x, REAL* y, int p);
        Interpolator(std::vector<REAL> &x, std::vector<REAL> &y);
        Interpolator(Interpolator& inter);
        ~Interpolator();
        void Clear();
        REAL Evaluate(REAL x);
        REAL Derivative(REAL x);    // 1st derivative
        REAL Derivative2(REAL x);   // 2nd derivative
        void SetMethod(INTERPOLATION_METHOD m);
        int Initialize();
        
        REAL MinX();
        REAL MaxX();

        REAL* GetXData();
        REAL* GetYData();
        gsl_spline* GetGslSpline();
        int GetNumOfPoints();
        INTERPOLATION_METHOD GetMethod();
        
        bool Freeze(); 
        void SetFreeze(bool f);
        void SetUnderflow(REAL min);
        void SetOverflow(REAL max);
        REAL UnderFlow();
        REAL OverFlow();

    private:
        INTERPOLATION_METHOD method;
        REAL* xdata, *ydata;
        bool allocated_data;    // true if we allocated xdata,ydata
        REAL minx,maxx;
        int points;
        bool ready;
        
        bool freeze;		// true if return freeze_under/overflow if
        REAL freeze_underflow;	// asked to evaluate interpolation
		REAL freeze_overflow;	// outside the spesified range
        
        // spline
        gsl_interp_accel *acc;
        gsl_spline *spline;

        // bspline
        gsl_bspline_workspace *bw;
        gsl_bspline_deriv_workspace *derbw;
        gsl_vector *B;
        gsl_vector *c;
        gsl_matrix *X;
        gsl_matrix *cov;
        gsl_multifit_linear_workspace *mw;

        static const int k=4;
        static const int ncoeffs = 12;
        static const int nbreak = ncoeffs-k+2;


};




#endif
