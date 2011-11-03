/*
 * Interpoate 2D data array
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _INTERPOLATION2D_H
#define _INTERPOLATION2D_H

/*
 * Interpolation helper
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
#include "config.hpp"
#include "tools.hpp"

#include "interpolation.hpp"

class Interpolator2D
{
    public:
        // data format [xind][yind]
        // Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
        Interpolator2D(std::vector<double>  &grid,
            std::vector< std::vector<double> > &data_);
        Interpolator2D(Interpolator2D& inter);
        ~Interpolator2D();
        void Clear();
        REAL Evaluate(REAL x, REAL y);
        REAL Derivative(REAL x, REAL y);    // 1st derivative
        REAL Derivative2(REAL x, REAL y);   // 2nd derivative
        void SetMethod(INTERPOLATION_METHOD m);
        int Initialize();

        INTERPOLATION_METHOD GetMethod();

    private:
        INTERPOLATION_METHOD method;
        REAL* ypoints;
        std::vector<Interpolator*> interpolators;

        bool ready;


};




#endif
