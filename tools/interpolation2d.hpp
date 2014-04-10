/*
 * Interpoate 2D data array
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _INTERPOLATION2D_H
#define _INTERPOLATION2D_H



#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <vector>
#include "config.hpp"
#include "tools.hpp"

#include "interpolation.hpp"

/**
 * 2D interpolator
 * Interpolates given data using spline (goes trough every data point)
 * or bspline (=noisy data)
 * Uses GSL
 *
 * Xdata and ydata pointers are saved, but not used for interpolation purposes
 * If user frees the allocated memory, one should be sure that these pointers
 * are not asked from this class!
 *
 * Uses Interpolator class to perform 1D interpolations
 */

class Interpolator2D
{
    public:
        // data format [xind][yind]
        // Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
        /**
         * Initialize 2D interpolator
         * Assumes that the interpolation grid is rectangular, and
         * coordinates are given in the grid vector.
         * Grid point corresponding to indexes (x,y) is
         * grid(x,y) = (grid[x], grid[y])
         * @param grid Vector of x and y coordinates
         * @param data 2D vector of datapoints, data[i][j] = f(grid[i], grid[j])
         * @see Interpolator
         */
        Interpolator2D(std::vector<double>  &grid,
            std::vector< std::vector<double> > &data_);
        Interpolator2D(Interpolator2D& inter);
        ~Interpolator2D();
        void Clear();
        double Evaluate(double x, double y);
        double Derivative(double x, double y);    // 1st derivative
        double Derivative2(double x, double y);   // 2nd derivative
        void SetMethod(INTERPOLATION_METHOD m);

        INTERPOLATION_METHOD GetMethod();

    private:
        INTERPOLATION_METHOD method;
        std::vector<double> ypoints;
        std::vector<Interpolator*> interpolators;

        bool ready;


};




#endif
