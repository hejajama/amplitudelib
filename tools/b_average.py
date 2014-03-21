#!/usr/bin/python
#-*- coding: UTF-8 -*-
# Computes impact parameter average
# \int \der^2 b f(b) / \int \der^2 b rho_inel(b)
# CLI arguments are fileprefix, filepostfix, bmin, bmax, bstep, protonb protonfile A model
#
# Filepostfix "-" means no postfix
#
# Part of AmplitudeLib 
# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013


# Configs
mbgevsqr = 2.568


import sys 
# filename is fileprefix+b+filepostfix
if (len(sys.argv)<10):
    print "Syntax: " + sys.argv[0] + " fileprefix filepostfix minb maxb bstep protonb protonfile A model"
    sys.exit(0)

fileprefix=sys.argv[1]
filepostfix=sys.argv[2]
if filepostfix=="-":
    filepostfix=""
minb=int(sys.argv[3])
maxb=int(sys.argv[4])
bstep=int(sys.argv[5])

switch_proton_b=int(sys.argv[6])
proton_file=sys.argv[7]

A=int(sys.argv[8])
sigmann=0
sqrts=int(sys.argv[10])
if sqrts==5020 or sqrts==7000:
    sigmann=70 * mbgevsqr
elif sqrts==8800:
    sigmann=75*mbgevsqr
elif sqrts==200:
    sigmann=42*mbgevsqr
else:
    sys.stderr.write("ERROR! Unknown A " + str(A) +"\n")
    sys.exit(-1)
    

model=sys.argv[9]
sigma02=0
if model=="mv1":
    sigma02=18.81*mbgevsqr
elif model=="mvgamma":
    sigma02=16.45*mbgevsqr
elif model=="mve":
    sigma02=16.36*mbgevsqr
elif model=="mvgamma_oma":
    sigma02=41.985
else:
    sys.stderr.write("ERROR! Unknown model " + model + "\n")
    

print "# b-average, minb " + str(minb) + " maxb " + str(maxb) + " bstep " + str(bstep) +" sigmann " + str(sigmann/mbgevsqr) +" mb, sigma02 " + str(sigma02/mbgevsqr)+ " mb, A " + str(A) +", switching to proton file " + proton_file + " at b>" + str(switch_proton_b) +" GeV^-1"


import os

sys.path.append("/nashome2/hejajama/lib/")
sys.path.append("/home/hejajama/lib/")
sys.path.append("../lib/")
import math
from matplotlibhelper import *
import pylab
import scipy.integrate
from scipy import interpolate
import locale
locale.setlocale(locale.LC_ALL,"C")

def RA(A):
    return fmgev*(1.12*pow(A, 1.0/3.0) - 0.86*pow(A, -1.0/3.0))

#def WoodsSaxon(r,A):
#     norm = 3.0/(4.0 * pi * pow(RA(A), 3)) * 1.0/(1.0 + pow( pi * ws_delta/RA(A),2))
#     return norm / (1.0 + exp((r-RA(A))/ws_delta) )
     
#def WSint(z, b, A):
#    return WoodsSaxon(sqrt(b*b+z*z), A)

#def TA(r, A):
#    result = scipy.integrate.quad(WSint, -999, 999, args=(r,A))
#    return result[0]

#def Inthelper_bint(b, interpolator):
#    return 2.0*pi*b*interpolabor(b)

def Nbin(b, A):
    return A*TA(b,A)*sigmann

# interpolators for every b
interpolators=[]

xvals=[]
bvals=[]
b=minb
while b <= min(switch_proton_b, maxb):
    fname = fileprefix + str(b) + filepostfix
    
    
    xdata=[]
    ydata=[]
        
    print "#" + fname
    try:
        readfile_xy(fname, xdata, ydata,xcol=0, ycol=1)
    except IOError:
        sys.stderr.write("ERROR: Nukefile " + fname + " does not exist\n")
        sys.exit(-1)
    except ValueError:
        sys.stderr.write("ERROR: Nukefile " + fname + " is invalid\n")
    #xdata=list(reversed(xdata))
    #ydata=list(reversed(ydata))
    if (b==minb):
        xvals=xdata
    try:
        interp=interpolate.interp1d(xdata, ydata, kind="cubic")
    except:
        sys.stderr.write("Interpolation failed, too few datapoints in file " + fname + "\n")
        sys.exit(0)
    interpolators.append(interp) 
    bvals.append(b)   
    b=b+bstep

# read proton file
proton_interpolator=None   
if float(switch_proton_b) < float(maxb):
    xdata=[]
    ydata=[]
    try:
        readfile_xy(proton_file, xdata, ydata, xcol=0, ycol=1)
    except IOError:
        sys.stderr.write("ERROR: Proton file " + proton_file + " does not exist\n")
        sys.exit(-1)
    except ValueError:
        sys.stderr.write("ERROR: Proton file " + proton_file +" is invalid\n")
        sys.exit(-1)
    
    try:
        proton_interpolator = interpolate.interp1d(xdata, ydata, kind="cubic")
    except ValueError:
        sys.stederr.write("ERROR: Too few proton datapoints (" + str(len(ydata)) +"), file " + proton_filename)
        sys.exit(-1)

# compute normalization \int der^2 b rho_inel(b)
normlist = []
normlist_b = []
b=minb
while b<=maxb:
    normlist_b.append(b)
    normlist.append(2.0*pi*b*rho_inel(b, A, sigmann))
    b=b+0.1
normalization = scipy.integrate.simps(normlist, normlist_b)

#nbin integral for the proton part
nbinvals=[]
nbinvals_blist=[]
b=switch_proton_b
while b<=maxb:
    nbinvals_blist.append(b)
    nbinvals.append(2.0*pi*b*Nbin(b, A))
    b=b+0.1
if len(nbinvals)==0:
    nbin_proton_coef=0
else:
    nbin_proton_coef = scipy.integrate.simps(nbinvals, nbinvals_blist)
print "# proton nbin coef: " + str(nbin_proton_coef)

for x in xvals:
    intvals=[]
    i=0
    for i in range(len(bvals)):
        try:
            intvals.append(interpolators[i](x)*2.0*pi*bvals[i])
        except ValueError as detail:
            sys.stderr.write("ValueError: x=" + str(x) + "\n")
            sys.exit(1)
        if (interpolators[i](x)*2.0*pi*bvals[i] < 0):
            print "#WWWTTTTTTTFFFFFFFFFF, b=" + str(bvals[i]) +" pt " + str(x) + " interp " + str(interpolators[i](x)*2.0*pi*bvals[i])
    
    #print bvals
    #print intvals
    bintres = scipy.integrate.simps(intvals, bvals)
    ppres=0
    if float(switch_proton_b) < float(maxb):
        try:
            ppres = sigma02/sigmann*nbin_proton_coef * proton_interpolator(x)
        except ValueError:
            sys.stderr.write("ValueError: x=" + str(x) + "\n")
            sys.exit(1)
    print "# pA contrib " + str(bintres) + " pp contrib " + str(ppres) +" normalization " + str(normalization)
    res = (bintres + ppres)/normalization
    print str(x) + " " + str(res)



