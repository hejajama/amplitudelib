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
    cerr << "Interpolator2D is already intialized when its construction is called! "
        << LINEINFO << endl;
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
            << "it in that interval! " << LINEINFO << endl;
        if (x<ypoints[0]) x=ypoints[0]*1.00001;
        if (y<ypoints[0]) y=ypoints[0]*1.00001;
        if (x>maxy) y=maxy*0.999999;
        if (y>maxy) x=maxy*0.999999;
    }
    
    REAL res; int status;

    // Construct one new interpolator, we already have interpolators[]
    // evaluated at each y in the grid
    std::vector<double> tmpdata;
    for (uint yind=0; yind<interpolators.size(); yind++)
    {
        tmpdata.push_back(interpolators[yind]->Evaluate(x));
    }
    Interpolator inter(ypoints, tmpdata);
    inter.SetMethod(INTERPOLATE_SPLINE);
    status = inter.Initialize();

    if (status)
    {
        cerr << "Interpolator initialization failed at " << LINEINFO << endl;
        return 0;
    }

    res = inter.Evaluate(y);


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

    if (grid.size()<4 or data.size()<4)
    {
        cerr << "Grid size " << grid.size() << " or data size " << data.size()
            << " are too small for interpolation! " << LINEINFO << endl;
        return;
    }
    
    // Construct interpolators in x direction, we get # of ypoints
    // interpolators which can be used to evaluate any point (x,y)
    ypoints.clear();
    for (uint i=0; i < grid.size(); i++)
    {
        ypoints.push_back(grid[i]);
        if (i>0)
        {
            if (grid[i-1]>=grid[i])
            {
                cerr << "Grid[" << i-1 << "]=" << grid[i-1]
                    << " >= Grid[" << i << "]=" << grid[i]
                    << ", gridsize=" << grid.size()
                    << ", can't initialize 2D interpolator! " << LINEINFO << endl;
                return;
            }
        }
    }
    std::vector<double> tmpdata;
    for (uint yind=0; yind < grid.size(); yind++)
    {
        tmpdata.clear();
        for (uint xind=0; xind<grid.size(); xind++)
        {
            tmpdata.push_back(data[xind][yind]);
        }
        Interpolator* interp = new Interpolator(grid, tmpdata);
        interpolators.push_back(interp);
        interp->SetMethod(INTERPOLATE_SPLINE);
        interp->Initialize();
    }
    

    ready=true;
}

void Interpolator2D::SetMethod(INTERPOLATION_METHOD m)
{
    method = m;
}

void Interpolator2D::Clear()
{
    ypoints.clear();
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
