/*
 * General purpose interpolation class
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011-2013
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
#include <vector>
#include <iostream>

using std::cout;
using std::endl;
using std::cerr;



enum INTERPOLATION_METHOD {
    INTERPOLATE_SPLINE,
    INTERPOLATE_BSPLINE
};

class Interpolator
{
    public:
        Interpolator(double* x, double* y, int p);
        Interpolator(std::vector<double> &x, std::vector<double> &y);
        Interpolator(Interpolator& inter);
        ~Interpolator();
        void Clear();
        double Evaluate(double x);
        double Derivative(double x);    // 1st derivative
        double Derivative2(double x);   // 2nd derivative
        void SetMethod(INTERPOLATION_METHOD m);
        int Initialize();
        
        double MinX();
        double MaxX();

        double* GetXData();
        double* GetYData();
        gsl_spline* GetGslSpline();
        int GetNumOfPoints();
        INTERPOLATION_METHOD GetMethod();
        
        bool Freeze(); 
        void SetFreeze(bool f);
        void SetUnderflow(double min);
        void SetOverflow(double max);
        double UnderFlow();
        double OverFlow();

    private:
        INTERPOLATION_METHOD method;
        double* xdata, *ydata;
        bool allocated_data;    // true if we allocated xdata,ydata
        double minx,maxx;
        int points;
        bool ready;
        
        bool freeze;		// true if return freeze_under/overflow if
        double freeze_underflow;	// asked to evaluate interpolation
		double freeze_overflow;	// outside the spesified range
        
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
