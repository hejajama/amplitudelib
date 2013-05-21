#!/usr/bin/python
#-*- coding: UTF-8 -*-
# Computes impact parameter average
# \int \der^2 b f(b) / \int \der^2 b rho_inel(b)
# CLI arguments are bmin and bmax
#
# Part of AmplitudeLib 
# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013


# Configs
fileprefix="sinc/lhc_pA/theory/paperi/minft0/pA_b"
filepostfix="_hybrid_y_2_sqrts_5020_cteq_dss_lo"
# filename is fileprefix+b+filepostfix

mbgevsqr = 2.568
sigmann = 60 * mbgevsqr
A=208

import sys 
minb=float(sys.argv[1])
maxb=float(sys.argv[2])
bstep=2



import os

sys.path.append("/nashome2/hejajama/lib/")
import math
from matplotlibhelper import *
import pylab
import scipy.integrate
from scipy import interpolate

def Inthelper_bint(b, interpolator):
    return 2.0*pi*b*interpolabor(b)

# interpolators for every b
interpolators=[]

xvals=[]
bvals=[]
b=minb
while b <= maxb:
    fname = fileprefix + str(b) + filepostfix
    print "#" + fname
    
    xdata=[]
    ydata=[]
    readfile_xy(fname, xdata, ydata)
    if (b==minb):
        xvals=xdata
    interp=interpolate.interp1d(xdata, ydata, kind="cubic")
    interpolators.append(interp) 
    bvals.append(b)   
    b=b+bstep


# compute normalization \int der^2 b rho_inel(b)
normlist = []
for b in bvals:
    normlist.append(2.0*pi*b*rho_inel(b, A, sigmann))

normalization = scipy.integrate.simps(normlist, bvals)

for x in xvals:
    intvals=[]
    i=0
    for i in range(len(bvals)):
        intvals.append(interpolators[i](x)[0]*2.0*pi*bvals[i])
    #print bvals
    #print intvals
    res = scipy.integrate.simps(intvals, bvals) / normalization
    print str(x) + " " + str(res)



