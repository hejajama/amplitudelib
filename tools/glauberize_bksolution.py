#!/usr/bin/python
#-*- coding: UTF-8 -*-

# Takes solution of the bk equation and glauberizes it with given b and A
# Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2013

import scipy
from math import*
import sys

scipy.pkgload()

A=208
b=0
b=float(sys.argv[2])


fmgev=5.0677
gamma=1.135
lambdaqcd=0.241
#sigma0=32.895*2.568 #mvgamma
sigma0=28.7977*2.568 #1507.03651
ws_delta = 0.54*fmgev

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



    
def GlauberN(n,b,A):
    return 1.0 - exp(-A*TA(b,A)*sigma0/2.0 * n )



# Read bk solution
params=[0,0,0,0]  # rmin  dr N x0
yvals=[]
nvals=[]  # list of lists, every list contains nvals amplitude values
f=open(sys.argv[1])
lines = f.readlines()
f.close()

parami=0
tmpnvals=[]
for l in lines:
    if l[:3]=="###":
        if (parami<4):
            params[parami]=float(l[3:])
            parami=parami+1  
        else:
            yvals.append(float(l[3:]))
            if (len(tmpnvals)>0):
                nvals.append(tmpnvals)
                tmpnvals=[]
    elif l[0]=="#":
        continue
    else:
        tmpnvals.append(float(l))
nvals.append(tmpnvals)

# Glauberize and print
print "# Glauberized from file " + sys.argv[1] + ", A=" + str(A) + ", b=" + str(b) +" sigma0 = " + str(sigma0) + " (in GeV^n)"
for p in params:
    print "###" + str(p)
for i in range(len(yvals)-1):
    print "###" + str(yvals[i])
    for n in nvals[i]:
        print GlauberN(n, b, A)





#print "#b = " + str(sys.argv[1])
#b=float(sys.argv[1])
#r=pow(10.0,-6)
#while(r<100):
#    print str(r) + " " + str(GlauberN(r,b,208))
#    r=r*1.0472494

