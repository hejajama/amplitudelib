#!/usr/bin/python
#-*- coding: UTF-8 -*-
# Create dipole amplitude for the nucleus at the initial condition
# using the procedure described in 
# T. Lappi, H. Mäntysaari, arXiv:1309.6963 [hep-ph]
# H. Mäntysaari <heikki.mantysaari@jyu.fi>, 2013-2014

import scipy
from math import*
import sys

scipy.pkgload()

fmgev=5.0677
lambdaqcd=0.241

ws_delta = 0.54*fmgev
nukeA=197

#mvgamma
#qsqr=0.165
#gamma=1.135
#ec=1
#sigma0=84.47436
#qsqr=0.1586
#gamma=1.1290
#ec=1
#sigma0=82.97

#mv
#qsqr=0.104
#gamma=1.0
#ec=1
#sigma0=96.586       # sigma_0 in [GeV^{-2}]

#mve
qsqr=0.060
ec=18.9
gamma=1
#sigma0=83.985
sigma0 = 16.36*2.0*2.568  # = 84.025

def RA(A):
    return fmgev*(1.12*pow(A, 1.0/3.0) - 0.86*pow(A, -1.0/3.0))

def WoodsSaxon(r,A):
     norm = 3.0/(4.0 * pi * pow(RA(A), 3)) * 1.0/(1.0 + pow( pi * ws_delta/RA(A),2))
     return norm / (1.0 + exp((r-RA(A))/ws_delta) )
     
def WSint(z, b, A):
    return WoodsSaxon(sqrt(b*b+z*z), A)

def TA(r, A):
    result = scipy.integrate.quad(WSint, -999, 999, args=(r,A))
    return result[0]


# Expanded dipole amplitude
def Np(r):
	return pow(r*r*qsqr,gamma)/4.0*log(1.0/(r*lambdaqcd)+ec*e)
	
# Dipole amplitude as a function of b, Eq. (17) from the Ref.
def GlauberN(r,b,A):
    return 1.0 - exp(-A*TA(b,A)*sigma0/2.0 * Np(r) )

print "#b = " + str(sys.argv[1]) + " sigma0 " + str(sigma0) + " gamma " + str(gamma) + " qsp^2 " + str(qsqr) + " ec " + str(ec) + " A " + str(nukeA)
print "# All dimensionfull units are in GeV^n"
b=float(sys.argv[1])
r=pow(10.0,-6)
while(r<100):
    print str(r) + " " + str(GlauberN(r,b,nukeA))
    r=r*1.018610170155976

