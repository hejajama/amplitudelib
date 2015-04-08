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
import pdb

def format_number(i, decimals=4):
    #return '%g' % (Decimal(str(i)))
    try:
        return ("%0." + str(decimals)+"f") % float(i)
    except:
        pdb.set_trace()

fmgev = 5.067

heraparams = 4

minsigma0 = 4.0*math.pi*2.0*3.5
maxsigma0 = 4.0*math.pi*2.0*4.5

minsigma0 = 40
maxsigma0 = 150

def PrintMatrix(m):
    for y in range(len(m)):
        row=""
        for x in range(len(m[y])):
            row = row + " " + str("%0.4f" % round(m[y][x],4))
        print row


# Determines the best normalization coefficient for the given file
# Returns a list where elements are 
# [coefficient, xsisqr]
# Coefficient is in GeV^-2

def FitHera(filename, coefstep=0.1, cc=-1):
    expdata=[]
    experr=[]
    theory=[]
    
    
    params=filename.split("_")
    #if (params[6]!="0.3152"):   # "correct" C^2
    #    return [0,999999999]
    #if (params[4]!="1" and params[8]!="1.dat"):
    #    return [0,999999999]  # \gamma != 1 AND ec != 1, dont like this...
    #if (params[8]!="1.dat"):
    #    return [0,999999]
    f = open(filename,"r")
    lines=f.readlines()
    n=len(lines)
    if (n<200):
        print "Skipping file " + filename + " as there are only " + str(n) + " datapoints!!!"
        return [0, 999999999]
    for i in range(n):
	    if (lines[i][0]=="#"):
	        continue
	    if (len(lines[i])<5):
	        continue
	    s=lines[i].split()
	    #if (len(s)<6):
	    #    print filename
	    #    print s
	    if (float(s[0])>50):
	        continue
	    expdata.append(float(s[3]))
	    experr.append(float(s[4]))
	    theory.append(float(s[5]))
    f.close()
    coef=minsigma0
    n=len(expdata)
    
    xsisqr=0
    bestfit_coef = coef
    bestfit_xsisqr = 999999999999;
    
    cc=float(cc)
    if (cc>0):
        coef=cc
        magsigma0=coef+0.000000000001
        coefstep=999
    
    while (coef <= maxsigma0):
        
        xsisqr = 0
        for i in range(len(expdata)):
            xsisqr = xsisqr + math.pow((expdata[i]-coef*theory[i])/experr[i],2)
        #xsisqr = xsisqr/float(n-heraparams)
        #if (xsisqr > 10*bestfit_xsisqr):
        #    coef=9999999
        if (xsisqr < bestfit_xsisqr):
            bestfit_xsisqr = xsisqr
            bestfit_coef = coef
        coef += coefstep
    
    #if (abs(1.0-bestfit_xsisqr) < 0.11):
        #print filename
    #print "file " + filename + ", \\sigma_0 = : " + str(bestfit_coef) + " GeV^-2 = " + str(round(bestfit_coef / (fmgev*fmgev) * 10.0,2)) + " mb, xsisqr = " + str(round(bestfit_xsisqr,3)) 
    
    xsisqr_dof = bestfit_xsisqr / float(len(expdata))
    
    # Print good fit result
    if (xsisqr_dof < 5):
        print filename + " xsisqr_dof " + str(xsisqr_dof) + " coef " + str(bestfit_coef) + " GeV^-2 "
    return [bestfit_coef, bestfit_xsisqr]

def ComputeHeraXsiSqr(fileprefix, params):
    qsqr=params[0]
    gamma=params[1]
    csqr=params[2]
    sigma = float(params[3])
    filename = fileprefix + "mv_qsqr_" + format_number(qsqr,2) + "_gamma_" + format_number(gamma,0) + "_csqr_" + format_number(csqr,0) + "_ec_1_freeze_0.dat"
    #print [filename,params]
    try:
        f = open(filename,"r")
    except:
        print "File " + filename + " does not exist"
        return 9999999
    lines=f.readlines()
    n=len(lines)
    expdata=[]
    experr=[]
    theory=[]
    if (n<200):
        print "Skipping file " + filename + " as there are only " + str(n) + " datapoints!!!"
        return [0, 999999999]
    for i in range(n):
	    if (lines[i][0]=="#"):
	        continue
	    if (len(lines[i])<5):
	        continue
	    s=lines[i].split()
	    if (float(s[0])>50):
	        continue
	    expdata.append(float(s[3]))
	    experr.append(float(s[4]))
	    theory.append(float(s[5]))
    f.close()
    xsisqr = 0
    for i in range(len(expdata)):
        try:
            xsisqr = xsisqr + math.pow((expdata[i]-sigma*theory[i])/experr[i],2)
        except:
            pdb.set_trace()
    #
    #xsisqr_dof = xsisqr/float(len(expdata)-heraparams)
    
    
    
    
    return xsisqr

# Find the best K factor to hadronic data   
# Returns [Kfactor, xsisqr] 
def FitHadron(filename, sigma0, coefstep=0.1):
    if not os.path.exists(filename):
        #print "Hadrondatafile " + filename + " does not exists! "
        return [0, 9999]


    expdata=[]
    experr=[]
    theory=[]
    
    sigma_inel_lhc = 60.0 * 0.1 * fmgev*fmgev

    f = open(filename,"r")
    lines=f.readlines()
    n=len(lines)
    for i in range(n):
        if (lines[i][0]=="#"):
            continue
        s=lines[i].split()
        expdata.append(float(s[4]))
        experr.append(float(s[5]))
        th = float(s[6])
        th = th * pow(sigma0/2.0,2) / sigma_inel_lhc
        theory.append(th)
    
    k=0.1
    bestfit_k=0
    bestfit_xsisqr=9999999
    while(k<20):
        xsisqr=0
        for i in range(len(expdata)):
            xsisqr = xsisqr + math.pow((expdata[i]-k*theory[i])/experr[i],2)
        xsisqr = xsisqr/float(n)
        #if (abs(1.0-xsisqr) > 2*abs(1.0-bestfit_xsisqr)):
        #    coef=9999999
        #if (abs(1.0-xsisqr) < abs(1.0-bestfit_xsisqr)):
        if (xsisqr < bestfit_xsisqr):
            bestfit_xsisqr = xsisqr
            bestfit_k = k
        k += coefstep
        
    return [bestfit_k, bestfit_xsisqr]
    
# Assume that there exists i s.t. l[i]==val
# Returns from the list l smallest possible value which is still larger than l[i]
def NextValue(l, val):
    mindiff=999
    mini = -1
    for i in range(len(l)):
        if (l[i] > val + 0.0000001):
            if (l[i]-val < mindiff):
                mini=i
                mindiff = l[i]-val
    return l[mini]   

def PrevValue(l, val):
    mindiff=999
    mini = -1
    for i in range(len(l)):
        if (l[i] < val - 0.00000001):
            if (val-l[i] < mindiff):
                mini=i
                mindiff = val-l[i]
    return l[mini]        





# Goes trough every datafile in folder fileprefix and computes chi^2
# Prints the best fit result and also other good fit results
# Filenames must be in the form
# text_qsqr_0.14_gamma_1_csqr_7_ec_1_freeze_0.dat
# where 0.14 is Qs^2, gamma=1, C^2=7 and e_c=1 
def CorrelationAnalysis():
    
    #Folder where files are read  
    fileprefix="mvfit/"
    
    # read data
    params=[]
    qsqrvals=[]
    gammavals=[]
    csqrvals=[]
    sigmavals=[]
    
    minxsisqr = 99999
    min_index=0
    sigmastep=0.01
    
    
    
    
    
    #fileprefix = "sinc/"
    for files in sorted( os.listdir(fileprefix) ):
        if files.endswith(".dat"):   
            fname=files.split("_")
            [sigma,xsisqr] = FitHera(fileprefix + files,sigmastep)
            #if (xsisqr < 1000):
            #    print files + " xsisqr " + str(xsisqr)  + " sigma " + str(sigma)
            if (xsisqr > 9999):
                continue
            
            # Only solutions with gamma=1, comment when analyzing aamqs
            if (fname[4] != "1"):
                continue
            qsqrvals.append(float(fname[2]))
            gammavals.append(float(fname[4]))
            csqrvals.append(float(fname[6]))
            
            sigmavals.append(sigma)
            if (xsisqr < minxsisqr):
                minxsisqr = xsisqr
                min_index = len(qsqrvals)-1
    
    tmpsigma=30
    while (tmpsigma<200):
        sigmavals.append(tmpsigma)
        tmpsigma=tmpsigma+sigmastep
    
    
    params.append(qsqrvals)
    params.append(gammavals)
    params.append(csqrvals)
    params.append(sigmavals)
    
    
    components=4
    tmpvec=[0]*components
    hess = [[0]*components for x in xrange(components)]
    tmpparams=[0]*components
    tmpparams[0]= qsqrvals[min_index]
    tmpparams[1] = gammavals[min_index]
    tmpparams[2] = csqrvals[min_index]
    tmpparams[3] = sigmavals[min_index]
    
    
    
    
    print "Bestfit: Q_s^2 " + str(qsqrvals[min_index]) + " c^2 " + str(csqrvals[min_index]) + " gamma " + str(gammavals[min_index]) + " sigma " + str(sigmavals[min_index]) + " GeV^-2 xsisqr " + str(ComputeHeraXsiSqr(fileprefix, tmpparams))
    
    
    return
    
    #### Correlation analyziz, not working correctly?????????????????????????
    for x in range(len(params)):
        for y in range(len(params)):
            tmpparamspp = list(tmpparams)
            tmpparamspm = list(tmpparams)
            tmpparamsmp = list(tmpparams)
            tmpparamsmm = list(tmpparams)

            
            tmpparamspp[x] = NextValue(params[x], tmpparams[x] ) #)params[x][min_index])
            tmpparamspp[y] = NextValue(params[y], tmpparams[y] ) #params[y][min_index])
            tmpparamspm[x] = NextValue(params[x], tmpparams[x] )# params[x][min_index])
            tmpparamspm[y] = PrevValue(params[y], tmpparams[y] )#params[y][min_index])
            tmpparamsmp[x] = PrevValue(params[x], tmpparams[x] )#params[x][min_index])
            tmpparamsmp[y] = NextValue(params[y], tmpparams[y] ) #params[y][min_index])
            tmpparamsmm[x] = PrevValue(params[x], tmpparams[x] )#params[x][min_index])
            tmpparamsmm[y] = PrevValue(params[y], tmpparams[y] )#params[y][min_index])
            
            if (abs( abs(NextValue(params[x], params[x][min_index]) - params[x][min_index]) - abs(params[x][min_index] - PrevValue(params[x], params[x][min_index]) ) ) > 0.0001):
                print "Non-uniform grid at component " + str(x) + " minvalue " + str(params[x][min_index]) + " +h " + str(NextValue(params[x], params[x][min_index]) ) + " -h " + str(PrevValue(params[x], params[x][min_index]))
                pdb.set_trace()
            

            if (x==y):
                fp = ComputeHeraXsiSqr(fileprefix, tmpparamspp) - minxsisqr
                f = 0
                fm = ComputeHeraXsiSqr(fileprefix, tmpparamsmm) - minxsisqr
                h = tmpparamspp[x] - params[x][min_index]
                if (abs( (tmpparamspp[x]-params[x][min_index]) - (params[x][min_index]-tmpparamsmm[x]) ) > 0.00001):
                    print "Non-uniform grid at component " + str(x) + " minvalue " + str(params[x][min_index]) + " +h " + str(NextValue(params[x], params[x][min_index]) ) + " -h " + str(PrevValue(params[x], params[x][min_index]))
                #print "plus: " + str(fp)  + " plusparams " + str(tmpparamspp) + " minus "  + str(fm) + " h " + str(h)
                hess[x][y] = 0.5* (fp + fm - 2.0*f) / (h*h)
                
                if (fp<0 or fm<0 or h<0):
                    pdb.set_trace()
                
                #print tmpparamspp,tmpparams
                #print ComputeHeraXsiSqr(fileprefix,tmpparamspp),ComputeHeraXsiSqr(fileprefix,tmpparams)
                #print tmpparamspp, tmpparamsmm,fp,fm
                print "# 1st derivative: " + str((fp - fm)/(2.0*h)) + " 2nd derivative " + str(hess[x][x])
                #pdb.set_trace()
            else:   
                print tmpparamspp, tmpparamspm, tmpparamsmp, tmpparamsmm
                fpp = ComputeHeraXsiSqr(fileprefix, tmpparamspp) - minxsisqr
                fpm = ComputeHeraXsiSqr(fileprefix, tmpparamspm) - minxsisqr
                fmp = ComputeHeraXsiSqr(fileprefix, tmpparamsmp) - minxsisqr
                fmm = ComputeHeraXsiSqr(fileprefix, tmpparamsmm) - minxsisqr

                
                h1 = tmpparamspp[x] - tmpparams[x]
                h2 = tmpparamspp[y] - tmpparams[y]
                
                #print x,y,h1,h2,fpp,fpm,fmp,fmm
                try:
                    hess[x][y] = 0.5 * (fpp - fpm - fmp + fmm) / (4.0 * h1 * h2)
                except:
                    pdb.set_trace()
                    print "wtf"
    
    # Write Hessian as a band matrix
    hessband = [[0]*components for x in xrange(components)]
    u=3
    for i in range(components):
        for j in range(components):
            if (i<=j):
                #pdb.set_trace()
                hessband[u+i-j][j] = hess[i][j]
    
    print "======================"
    print "Hessian"
    print "qsqr  gamma   csqr  sigma"  
    print "=========================="
    PrintMatrix(hess)
    
    print "Eigenvals and eigenvecs"
    print "=========================="
    eigenvals,eigenvecs = linalg.eig(hess)
#    evals2 = linalg.eigenvals_banded(hessband)
  #  print eigenvals
 #   print evals2
    
    basis=[]
    for i in range(len(eigenvals)):
        basis.append(eigenvecs[:,i])    
    
    for i in range(len(eigenvals)):
        print eigenvals[i]
        print basis[i]
    

    dzvecs_p=[[0]*components for x in xrange(components)]
    dzvecs_m=[[0]*components for x in xrange(components)]
            
    dzmultipliers=[1,2,3,4,5,6]
    
    print ""
    print "===="
    print "Original params"
    print tmpparams
    
    for deltaz in dzmultipliers: # [1,2,3,4,5,6]:
        for i in range(components):
            dzvecs_p[i][i] = 1.0*deltaz
            dzvecs_m[i][i] = -1.0*deltaz
        
        #dzvecs_m[3][3]=-0.5*deltaz
        
        #print dzvecs_p[0]
        #print dzvecs_m[0]
        #print ""
        #print dzvecs_p[1]
        #print dzvecs_m[1]
        
        errorparams_p = [[0]*components for x in xrange(components)]
        errorparams_m = [[0]*components for x in xrange(components)]
        for eset in range(components):
            for i in range(components):
                sp=0
                sm=0
                for j in range(components):
                    sp = sp + basis[j][i] * math.sqrt(1.0/eigenvals[j]) * dzvecs_p[eset][j]
                    sm = sm + basis[j][i] * math.sqrt(1.0/eigenvals[j]) * dzvecs_m[eset][j]
                
                errorparams_p[eset][i] = params[i][min_index] + sp
                errorparams_m[eset][i] = params[i][min_index] + sm
                #print [errorparams_p[set][i], errorparams_m[set][i] ]
        
        print "#Errorparams (deltaz=" + str(deltaz) + ")"
        print "#==================="
        for set in range(components):
            plist = [ '%.4f' % elem for elem in errorparams_p[set] ]
            tmpxsisqr=ComputeHeraXsiSqr(fileprefix + "/errorset/", plist)
            #print "# " + str(plist)
            print str(deltaz) + " " + str(tmpxsisqr-minxsisqr)
            #print plist, tmpxsisqr-minxsisqr
            print str(plist[0]) + " " + str(plist[1]) + " 0.01 " + str(plist[2]) + " 1 0 fit/aamqs/errorset/mv_qsqr_" + str(plist[0]) + "_gamma_" + str(plist[1]) +"_csqr_"+str(plist[2])+"_ec_1_freeze_0.dat " + str(plist[3])
            
            plist = [ '%.4f' % elem for elem in errorparams_m[set] ]
            tmpxsisqr=ComputeHeraXsiSqr(fileprefix + "/errorset/", plist)
            #print "# " + str(plist)
            print str(deltaz) + " " + str(tmpxsisqr-minxsisqr)
            #print  plist, tmpxsisqr-minxsisqr
            print str(plist[0]) + " " + str(plist[1]) + " 0.01 " + str(plist[2]) + " 1 0 fit/aamqs/errorset/mv_qsqr_" + str(plist[0]) + "_gamma_" + str(plist[1]) +"_csqr_"+str(plist[2])+"_ec_1_freeze_0.dat " + str(plist[3])
            
            #print errorparams_p[set]
            #paramlist = [ '%.4f' % elem for elem in errorparams_p[set] ]
            #print ComputeHeraXsiSqr(fileprefix + "/errorset/", paramlist)
            #print errorparams_m[set]
            #print ""
    


# If we got a filename as a cli parameter, fit it to HERA
if len(sys.argv)>1:
    [normalization, xsisqr] = FitHera(sys.argv[1])
    print "Normalization " + str(normalization) + " GeV^(-2) = " + str(round(normalization / (fmgev*fmgev) * 10.0,2)) + " mb, xsisqr = " + str(round(xsisqr,3)) 
    if (len(sys.argv)>2):
        [normalization, xsisqr] = FitHera(sys.argv[1],1,sys.argv[2])
        print "Normalization " + str(normalization) + " GeV^(-2) = " + str(round(normalization / (fmgev*fmgev) * 10.0,2)) + " mb, xsisqr = " + str(round(xsisqr,3)) 

    sys.exit(0)
    
    
CorrelationAnalysis()




