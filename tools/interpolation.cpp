/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation.hpp"


/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator::Initialize()
{
    int status=0;
    if (ready)
    {
        cerr << "Interpolator is already initialized, clear it before re-initializing! "
            << LINEINFO << endl;
        return -1;
    }
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            acc = gsl_interp_accel_alloc();
            spline = gsl_spline_alloc(gsl_interp_cspline, points);
            status = gsl_spline_init(spline, xdata, ydata, points);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_vector *x = gsl_vector_alloc(points);
            gsl_vector *y = gsl_vector_alloc(points);
            gsl_vector *w = gsl_vector_alloc(points);

            for (int i=0; i< points; i++)
            {
                gsl_vector_set(x, i, xdata[i]);
                gsl_vector_set(y, i, ydata[i]);
                gsl_vector_set(w, i, 1.0);
            }
     
            /* allocate a cubic bspline workspace (k = 4) */
            bw = gsl_bspline_alloc(k, nbreak);
            derbw = gsl_bspline_deriv_alloc(k);
            B = gsl_vector_alloc(ncoeffs);
       
            X = gsl_matrix_alloc(points, ncoeffs);
            c = gsl_vector_alloc(ncoeffs);
       
            cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
            mw = gsl_multifit_linear_alloc(points, ncoeffs);
     
     
            // use uniform breakpoints
            gsl_bspline_knots_uniform(xdata[0], xdata[points-1], bw);
     
            /* construct the fit matrix X */
            for (int i = 0; i < points; ++i)
            {
               double xi = gsl_vector_get(x, i);
             
               /* compute B_j(xi) for all j */
               gsl_bspline_eval(xi, B, bw);
             
               /* fill in row i of X */
               for (int j = 0; j < ncoeffs; ++j)
               {
                  double Bj = gsl_vector_get(B, j);
                  gsl_matrix_set(X, i, j, Bj);
               }
            }
     
            /* do the fit */
            REAL chisq;
            gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

            gsl_vector_free(x);
            gsl_vector_free(y);
            gsl_vector_free(w);

            break;
    }
    ready=true;
    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        return -1;
    }
    return 0;   //ok, there is no error handling at the moment...
}


REAL Interpolator::Evaluate(REAL x)
{
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator::Initialize()?" << endl;
        return 0;
    }

    if (x<minx or x>maxx)
    {
		if (x < 0.9999*minx or x > 1.00001*maxx)	// if not true, no need to display error
			cerr << "x=" << x << " is not within limits [" << minx << ", " << maxx << "], forcing "
				<< "it in that interval! " << LINEINFO << endl;
        if (x<minx) x=minx*1.00001;
        if (x>maxx) x=maxx*0.999999;
    }
    
    REAL res, yerr; int status;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_e(spline, x, acc, &res);
            if (status)
            {
                cerr << "Interpolation failed at " << LINEINFO << ", error " << gsl_strerror(status)
                 << " (" << status << "), x=" << x << ", minx=" << xdata[0]
                 << ", maxx=" << xdata[points-1] << ", result=" << res << endl;
                 exit(1);
            }
            return res;
            break;
        case INTERPOLATE_BSPLINE:
            gsl_bspline_eval(x, B, bw);
            gsl_multifit_linear_est(B, c, cov, &res, &yerr);

            /*if (std::abs(yerr/res)>0.05 )
            {
                cerr << "Interpolation failed at " << LINEINFO << ": bspline result "
                << res << " pm " << yerr << " relerr " << std::abs(yerr/res) << endl;
            }*/
            return res;
            break;
    }

    cerr << "Interpolation method is invalid! " << LINEINFO << endl;
    return 0;   //Shoudn't end up here
}

REAL Interpolator::Derivative(REAL x)
{
    REAL res=0; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_matrix* mat = gsl_matrix_alloc(nbreak+k-2, 2);
            gsl_bspline_deriv_eval(x, 1, mat, bw, derbw);
            for (int i=0; i<ncoeffs; i++)
            {
                res += gsl_vector_get(c, i)*gsl_matrix_get(mat, i, 1);
            }
            gsl_matrix_free(mat);
            return res;
    }
    if (status)
        cerr << "An error occurred while evaluating the derivative at x=" << x
        << " result " << res << " " << LINEINFO << endl;

    return res;
}

REAL Interpolator::Derivative2(REAL x)
{
    REAL res; int status=0;
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            status = gsl_spline_eval_deriv2_e(spline, x, acc, &res);
            break;
        case INTERPOLATE_BSPLINE:
            cerr << "2nd derivative is not implemented for BSPLINE interpolation!"
            << " " << LINEINFO << endl;
            break;
    }

    if (status)
    {
        cerr << "2nd derivative interpolation failed at x=" << x <<", result "
        << res << " " << LINEINFO << endl;
    }
    return res;

}

Interpolator::Interpolator(REAL *x, REAL *y, int p)
{
    points=p;
    xdata=x;
    ydata=y;
    minx=x[0];
    maxx=x[p-1];
    method = INTERPOLATE_SPLINE;
    allocated_data=false;
    ready=false;

    for (int i=0; i<p; i++)
    {
        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
            }
        }
    }
}

Interpolator::Interpolator(std::vector<REAL> &x, std::vector<REAL> &y)
{
    points = x.size();
    xdata = new REAL[points];
    ydata = new REAL[points];
    allocated_data=true;

    for (uint i=0; i<x.size(); i++)
    {
        xdata[i]=x[i];
        ydata[i]=y[i];

        // Check that x values are monotonically increasing
        if (i>0)
        {
            if (xdata[i-1]>=xdata[i])
            {
                cerr << "Grid points are not monotonically increasing! grid["
                    << i-1 <<"]=" << xdata[i-1] <<", grid["<<i<<"]="<< xdata[i]
                    << " " << LINEINFO << endl;
            }
        }
    }
    minx=xdata[0]; maxx=xdata[x.size()-1];
    method = INTERPOLATE_SPLINE;
    ready=false;
}

void Interpolator::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
}

void Interpolator::Clear()
{
    switch(method)
    {
        case INTERPOLATE_SPLINE:
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
            break;
        case INTERPOLATE_BSPLINE:
            gsl_bspline_free(bw);
            gsl_bspline_deriv_free(derbw);
            gsl_vector_free(B);
            gsl_matrix_free(X);
            gsl_vector_free(c);
            
            gsl_matrix_free(cov);
            gsl_multifit_linear_free(mw);
            break;
    }

    if (allocated_data)
    {
        delete[] xdata;
        delete[] ydata;
    }
}

Interpolator::~Interpolator()
{
    Clear();

}

REAL* Interpolator::GetXData()
{
    return xdata;
}
REAL* Interpolator::GetYData()
{
    return ydata;
}
int Interpolator::GetNumOfPoints()
{
    return points;
}
INTERPOLATION_METHOD Interpolator::GetMethod()
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator::Interpolator(Interpolator& inter)
{
    points=inter.GetNumOfPoints();
    xdata = new REAL[points];
    ydata = new REAL[points];
    allocated_data=true;


    gsl_spline *tmpspline = inter.GetGslSpline();
    for (int i=0; i<points; i++)
    {
        xdata[i] = tmpspline->x[i];
        ydata[i] = tmpspline->y[1];
    }
    minx = xdata[0]; maxx=xdata[points-1];
    method = inter.GetMethod();
    ready=false;
    Initialize();
}

gsl_spline* Interpolator::GetGslSpline()
{
    return spline;
}

REAL Interpolator::MinX()
{
	return minx;
}

REAL Interpolator::MaxX()
{
	return maxx;
}
