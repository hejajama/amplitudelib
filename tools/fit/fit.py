# utf-8
# Compute chi^2 by comparing with HERA \sigma_r data
# Usage:
# fit.py datafile,
# where datafile is computed using f2fit program.
# Returns values for \sigma_0 and chi^2


import os
import math
import sys
from decimal import Decimal
from numpy import linalg
from scipy.optimize import minimize_scalar
import pdb

def format_number(i, decimals=4, sign=True):
    try:
        ss="%"
        if sign:
            ss = ss + "+"
        ss = ss + "08." + str(decimals)+"f"
        return ss % float(i)
    except:
        pdb.set_trace()

fmgev = 5.067


############ Convert filenames to parameter arrays and vice versa
def ParamsToFilename(params):
    fname = "mv_qsqr_" + format_number(params[0]) + "_gamma_" + format_number(params[1]) + "_lncsqr_" \
                + format_number(params[2]) + "_lnec_" + format_number(params[3]) +"_freeze_" + str(params[4]) +".dat"
    #fname = "mv_qsqr_" + ("%04.4f") % float(params[0]) + "_gamma_" + ("%07.4f") % float(params[1]) + "_lncsqr_" \
    #            + format_number(params[2]) + "_lnec_" + format_number(params[3]) +"_freeze_" + str(params[4]) +".dat"
    return fname

def FilenameToParams(fname):
        'Convert filename to parameter array'
        name = fname.split("_")
        return [name[2], name[4], name[6], name[8], "0"]
#############################


class FitResult:
    def __init__(self):
        self.sigma0=0
        self.chisqr=0
        self.chsqr_n = 0    # chi^2/N
        self.n = 0  # number of data points
        self.filename = ""


class FitFile:
    'Class to handle a single set of parameters'
    'Data can be divided into multiple different files'
    'Filenames must be in the standard format, but they can be located in different directories'

    def __init__(self):
        self.expdata            = []    # experimental sigma_r in GeV^(-2)
        self.experr             = []    # experimental error
        self.theory_light       = []    # light quark contribution 
        self.theory_c           = []    # c contribution
        self.theory_b           = []    # b contribution
        self.datapoints         = 0     # Number of data points
        
        self.maxqsqr            = 50    # Q^2 limits
        self.minqsqr            = 0     

        self.minpoints          = 10    # Must be at least this many data points

        self.dirs               = []    # Directories
        self.quarks             = []    # [light, charm, beauty] for all dirs

    def AddDir(self, d, light=True, charm=True, bottom=True):
        self.dirs.append(d)
        self.quarks.append( [light, charm, bottom] )
    
    def ReadData(self, params):
        ' Read datafile, parameters array is in the form'
        ' [Q^2, gamma, ln C^2, ln ec, freeze]'
        ' Filename example: mv_qsqr_0.0600_gamma_01.0000_lncsqr_+02.0000_lnec_+02.0000_freeze_0.dat'

        for d,q in zip(self.dirs, self.quarks):
            fname = ParamsToFilename(params)
            self.ReadDataFile(d + "/" + fname, q[0], q[1], q[2])

    def ReadDataFile(self, fname, light=True, charm=True, bottom=True):
        lines=[]
            
        try:
            f = open(fname,"r")
            lines=f.readlines()
        except IOError:
            print "Could not open file " + fname + ", file does not exist?"
            sys.exit(-1)
            
        for i in range(len(lines)):
            if (lines[i][0]=="#"):
                continue
            if (len(lines[i])<5):
                continue
            s=lines[i].split()
            try:
                if (float(s[0])>self.maxqsqr or float(s[0])<self.minqsqr):
                    continue
            except:
                continue    # Probably the last line has some cpu usage/timing information?
                 
            self.expdata.append(float(s[3]))
            self.experr.append(float(s[4]))

            if light:
                self.theory_light.append(float(s[5]))
            else:
                self.theory_light.append(0)
            if charm:
                self.theory_c.append(float(s[6]))
            else:
                self.theory_c.append(0)
            if bottom:
                self.theory_b.append(float(s[7]))
            else:
                self.theory_b.append(0)
        
        f.close()
        n=len(self.expdata)
        self.datapoints = n
        if (n < self.minpoints):
            raise NameError("Too few data points, only " +str(n) +" in the file " + str(fname) + ", should be at least " + str(minpoints))

        
    def ChiSqr(self, sigma0):
        'Compute chi^2 for given sigma0 (assumed to be in GeV^(-2)'
        chisqr=0
        
        for i in range(len(self.expdata)):
            chisqr += math.pow( (self.expdata[i] 
                -sigma0*( self.theory_light[i] + self.theory_c[i] + self.theory_b[i] ) )
                / self.experr[i] , 2)

        return chisqr

    def MinimizeChiSqr(self):
        fit = minimize_scalar(self.ChiSqr)

        result = FitResult
        result.sigma0 = fit.x
        result.chisqr = fit.fun
        result.chisqr_n = fit.fun / len(self.expdata)
        result.n = len(self.expdata)
        return result



class FitDir:
    'Fit sigma0 for every file in the directory separately and find the best parametrization'

    
    def __init__(self):
        #Settings
        self.show_good_fits = False     # Print to stdout good fits
        self.good_fit_limit = 5         # limit for chi^2/N

        # Data arrays
        self.filedirs       =[]
        self.quarks         =[]
        self.qsqrvals       =[]
        self.gammavals      =[]
        self.lncsqrvals     =[]
        self.lnecvals       =[]
        self.sigmavals      =[]
        self.chisqrvals     =[]

    def AddDir(self, d, light=True, charm=True, bottom=True):
        self.filedirs.append(d)
        self.quarks.append([light, charm, bottom])

    def Fit(self):
        # Go trough every file
        for f in sorted( os.listdir(self.filedirs[0]) ):
            params = FilenameToParams(f)

            fitter = FitFile()

            for d,q in zip(self.filedirs, self.quarks):
                fitter.AddDir(d, q[0], q[1], q[2])
            
            fitter.ReadData(params)
            fitres = fitter.MinimizeChiSqr()

            self.qsqrvals.append(params[0])
            self.gammavals.append(params[1])
            self.lncsqrvals.append(params[2])
            self.lnecvals.append(params[3])
            self.sigmavals.append(fitres.sigma0)
            self.chisqrvals.append(fitres.chisqr)

            if self.show_good_fits and fitres.chisqr_n < self.good_fit_limit:
                par = params
                par.append(fitres.sigma0)
                print "Good fit with params " + str(par)     +": chi^2/N = " + str(fitres.chisqr_n)

            #print "Parameters " + str(params) +" fit result " + str(fitres.chisqr)

        # Done, find the smallest chisqr
        minchisqr=9999999999999
        min_index=-1
        for i in range(len(self.chisqrvals)):
            if self.chisqrvals[i] < minchisqr:
                minchisqr = self.chisqrvals[i]
                min_index = i

        res = FitResult
        res.sigma0 = self.sigmavals[min_index]
        res.chisqr = self.chisqrvals[min_index]
        res.filename = ParamsToFilename([self.qsqrvals[min_index], self.gammavals[min_index],
            self.lncsqrvals[min_index], self.lnecvals[min_index],0])
        return res

    def DirInfo(self):
        'Return string containing list of dirs and settings'
        ss=""
        for d,q in zip(self.filedirs, self.quarks):
            ss=ss + "dir: " + d +", quarks: " + str(q) + "\n"
        return ss
        




if __name__ == "__main__":

    if len(sys.argv)==1:
        sys.argv.append("-help")
        
    if sys.argv[1]=="-help":
        print "HELP:"
        print "-file fname true/false true/false true/false: fit given file using quarks light/c/b"
        print "-dir path true/false true/false true/false: add path to global fit"

    if sys.argv[1]=="-file":
        # Fit single file
        light=True
        charm=True
        bottom=True
        if sys.argv[3]=="false":
            light=False
        elif sys.argv[3]!="true":
            print "Unknown argument " + sys.argv[3]
            sys.exit(-1)
        if sys.argv[4]=="false":
            charm=False
        elif sys.argv[4]!="true":
            print "Unknown argument " + sys.argv[4]
            sys.exit(-1)
        if sys.argv[5]=="false":
            bottom=False
        elif sys.argv[5]!="true":
            print "Unknown argument " + sys.argv[5]
            sys.exit(-1)

        fname=sys.argv[2]

        print "# Fitting file " + fname
        print "# Light quarks: " + str(light)
        print "# Charm: " + str(charm)
        print "# Bottom: " + str(bottom)

        fit=FitFile()
        fit.ReadDataFile(fname, light, charm, bottom)
        bestfit=fit.MinimizeChiSqr()
        print "#Best sigma0: " + str(bestfit.sigma0) +", chi^2=" + str(bestfit.chisqr) +", chi^2/N = " + str(bestfit.chisqr_n) +", N = " + str(bestfit.n)
        print bestfit.chisqr_n

        sys.exit(0)

    # add multiple paths to fit
    fitter = FitDir()
    verbose=False
    for i in range(len(sys.argv)):
        if sys.argv[i]=="-dir":
            light=True
            charm=True
            bottom=True
            if sys.argv[i+2]=="false":
                light=False
            elif sys.argv[i+2]!="true":
                print "Unknown argument " + sys.argv[i+2]
                sys.exit(-1)
            if sys.argv[i+3]=="false":
                charm=False
            elif sys.argv[i+3]!="true":
                print "Unknown argument " + sys.argv[i+3]
                sys.exit(-1)
            if sys.argv[i+4]=="false":
                bottom=False
            elif sys.argv[i+4]!="true":
                print "Unknown argument " + sys.argv[i+4]
                sys.exit(-1)

            fitter.AddDir(sys.argv[i+1], light, charm, bottom)
        if sys.argv[i]=="-v":
                verbose=True

    print "Running fit"
    print fitter.DirInfo()

    if verbose:
        fitter.show_good_fits = True
        fitter.good_fit_limit = 4

    res = fitter.Fit()
    print "Bestfit: " + res.filename + ", sigma0=" + str(res.sigma0) + ", chi^2=" + str(res.chisqr) 

    #fit = FitFile()
    #fit.AddDir("./data/")
    #fit.ReadData([0.06,1,2,2,0])



    #bestfit = fit.MinimizeChiSqr()
    #print "Best sigma0: " + str(bestfit.sigma0) +", chi^2=" + str(bestfit.chisqr) +", chi^2/N = " + str(bestfit.chisqr_n)

    #fitter = FitDir()
    #fitter.AddDir("/space/hejajama/hera/2015/data/mve/charmmassfit/mc_1.27/sigmar/")
    #fitter.AddDir("/space/hejajama/hera/2015/data/mve/charmmassfit/mc_1.27/sigmar_cc/", light=False, charm=True, bottom=False)
    #fitter.AddDir("/space/hejajama/hera/herafit_systemaattisesti/sigmar/mve/", light=True, charm=False, bottom=False)
    #res = fitter.Fit()
    #print "Bestfit: " + res.filename + ", sigma0=" + str(res.sigma0) + ", chi^2=" + str(res.chisqr) 


