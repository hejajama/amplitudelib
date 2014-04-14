#!/usr/bin/python
# Create BK solution file for the MV model with x-dependent Qs:
# Q_s = Sqrt(0.2) * (x0/x)^0.1, where x0=0.1
# N(r,x) = 1-exp[ -(r * Q_s)^2/4 * Log( 1/(r*0.241) + E ) )
# Prints the BK file to standrard output

import math

x0=0.01
minr=1e-6
rpoints=400
rstep=1.047249412298149
maxy=16
ystep=0.2

def Qs(x):
    return math.sqrt(0.2) * math.pow(x0/x, 0.1)

def N(r, x):
    return 1.0 - math.exp( - math.pow(r*Qs(x), 2)/4.0 * math.log( 1/(r*0.241) + math.e) )

print "# Dipole amplitude created from standard MV model with Q_s = Sqrt(0.2) * (x0/x)^0.1"
print "# N(r,x) = 1-exp[ -(r * Q_s)^2/4 * Log( 1/(r*0.241) + E ) )"
print "#"
print "###" + str(minr)
print "###" + str(rstep)
print "###" + str(rpoints)
print "###" + str(x0)

y=0
while(y <=16.0):
    print "###" + str(y)
    x = x0*math.exp(-y)
    for rind in range(rpoints):
        r = minr*math.pow(rstep, rind)
        print N(r, x)
    y+=ystep


