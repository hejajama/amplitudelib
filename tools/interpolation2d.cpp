/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "interpolation2d.hpp"


/*
 * Intialize interpolation
 * Returns -1 in case of error, 0 otherwise
 */
int Interpolator2D::Initialize()
{
    if (ready)
    {
        cerr << "Can't initialize Interpolator2D, already initialized?" << LINEINFO<< endl;
        return -1;
    }
    
    for (uint i=0; i<interpolators.size(); i++)
    {
        interpolators[i]->SetMethod(method);
        int status = interpolators[i]->Initialize();
        if (status)
        {
            cerr << "Interpolator initialization failed at " << LINEINFO << endl;
            Clear();
            return -1;
        }
    }

    ready=true;
    return 0;
}


REAL Interpolator2D::Evaluate(REAL x, REAL y)
{
    if (!ready)
    {
        cerr << "Interpolator is not ready! Did you forget to call Interpolator2D::Initialize()?" << endl;
        return 0;
    }

    REAL maxy = ypoints[interpolators.size()-1];

    if (x<ypoints[0] or x>maxy or y<ypoints[0]
        or y>maxy)
    {
        cerr << "x or y is not within limits [" << ypoints[0] << ", " <<
            maxy << "], forcing "
            << "it in that interval!" << endl;
        if (x<ypoints[0]) x=ypoints[0]*1.00001;
        if (y<ypoints[0]) y=ypoints[0]*1.00001;
        if (x>maxy) y=maxy*0.999999;
        if (y>maxy) x=maxy*0.999999;
    }
    
    REAL res; int status;

    // Construct one new interpolator, we already have interpolators[]
    // evaluated at each y in the grid
    REAL* tmpdata = new REAL[interpolators.size()];
    for (uint yind=0; yind<interpolators.size(); yind++)
    {
        tmpdata[yind] = interpolators[yind]->Evaluate(x);
    }
    Interpolator inter(ypoints, tmpdata, interpolators.size());
    inter.SetMethod(INTERPOLATE_SPLINE);
    status = inter.Initialize();

    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        delete[] tmpdata;
        return 0;
    }

    res = inter.Evaluate(y);

    delete[] tmpdata;

    return res;
}

REAL Interpolator2D::Derivative(REAL x, REAL y)
{
    cerr << "2D interpolator derivative is not implemented! " << LINEINFO << endl;
    return 0;
}

REAL Interpolator2D::Derivative2(REAL x, REAL y)
{
    cerr << "2D interpolator derivative is not implemented! " << LINEINFO << endl;
    
    return 0;

}

// Initialize, data format: [xind][yind]
// Grid is rectangular, that is, grid(x,y) = (grid[x], grid[y])
Interpolator2D::Interpolator2D(std::vector<double>  &grid,
            std::vector< std::vector<double> > &data)
{
    // Construct interpolators in x direction, we get # of ypoints
    // interpolators which can be used to evaluate any point (x,y)

    ypoints = new double[grid.size()];
    for (uint i=0; i < grid.size(); i++)
    {
        ypoints[i] = grid[i];
        if (i>0)
        {
            if (ypoints[i-1]>=ypoints[i])
            {
                cerr << "Grid[" << i-1 << "]=" << ypoints[i-1]
                    << " >= Grid[" << i << "]=" << ypoints[i]
                    << ", gridsize=" << grid.size()
                    << ", can't initialize 2D interpolator! " << LINEINFO << endl;
                return;
            }
        }
    }
    double *tmpdata = new double[grid.size()];
    for (uint yind=0; yind < grid.size(); yind++)
    {
        for (uint xind=0; xind<grid.size(); xind++)
        {
            tmpdata[xind] = data[xind][yind];
        }
        Interpolator* interp = new Interpolator(ypoints, tmpdata, grid.size());
        interpolators.push_back(interp);
        interp->SetMethod(INTERPOLATE_SPLINE);
        interp->Initialize();
    }
    
    delete[] tmpdata;

    ready=true;
}

void Interpolator2D::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
}

void Interpolator2D::Clear()
{
    delete[] ypoints;
    for (uint i=0; i<interpolators.size(); i++)
    {
        delete interpolators[i];
    }
    interpolators.clear();

    ready=false;
}

Interpolator2D::~Interpolator2D()
{
    Clear();

}


INTERPOLATION_METHOD Interpolator2D::GetMethod()
{
    return method;
}

// Copy data from given class and initialize this, as this is
// the copy constructor
Interpolator2D::Interpolator2D(Interpolator2D& inter)
{
    cerr << "Interpolator2D copy constructor may not work???? " << LINEINFO << endl;
}
