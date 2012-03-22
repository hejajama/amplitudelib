/*
 * BK equation solver
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2011
 */

#ifndef _DATAFILE_H
#define _DATAFILE_H

#include "../tools/tools.hpp"
#include "../tools/config.hpp"
#include <sstream>
#include <vector>
/* Read data from datafiles, file format is defined in file README
 */

class DataFile
{
    public:
        DataFile(string fname);
        double MinR();
        double RMultiplier();
        int RPoints();
        double MaxY();
        double X0();

        void GetData(std::vector< std::vector<double> > &n,
            std::vector<double> &rapidities);

    private:
        string filename;
        std::vector<std::vector <double> > data;
        std::vector<double> yvals;
        double minr;
        double r_multiplier;
        int rpoints;
        double x0;

};

#endif
