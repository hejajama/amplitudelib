/*
 * Interpolator
 * Small example program/tool which uses the Interpolator class
 * Reads data from standard input and prints a finer grid to the
 * standard output generated using the SPLINE/BSPLINE interpolation
 * implemented in the Interpolator class
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2012
 */

#include "interpolation.hpp"
#include "config.hpp"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

int main(int argc, char* argv[])
{
    // Read all standard input to tables
    vector<double> xdata, ydata;
    string line;
    while (cin){
        getline(cin, line);
        if (line.size()<3 or line[0]=='#')
            continue;
        stringstream ss(line);
        double x,y;
        ss >> x; ss>>y;
        xdata.push_back(x); ydata.push_back(y);
    }

    // Generate interpolator
    Interpolator interp(xdata, ydata);
    interp.SetMethod(INTERPOLATE_SPLINE);
    interp.Initialize();

    // Print grid fined 10x
    double minx=xdata[0];
    double maxx=xdata[xdata.size()-1];
    int points = xdata.size()*10;

    for (int i=0; i < points; i++)
    {
        double tmpx = minx + (double)i/((double)points-1.0)*(maxx-minx);
        cout << tmpx << " " << interp.Evaluate(tmpx) << endl;
    }

    return 0;

}
