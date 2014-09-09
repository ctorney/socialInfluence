
#!/usr/bin/python

import sympy as sp
import scipy as sc
from scipy import optimize
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import math as m
import matplotlib as mpl

K = 80       
wg = 0.3

ws = 0.514
alpha = 0.0
NA = 64.0
N = 1024
def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*m.erf(m.sqrt(wg)*(1.0-gc))
    

def tup( xx ):
    return psw(K*xx)

def tdown( x,xx ):
    return (x) * (1 -  psw(K*xx))

numX = 64
Xs = np.zeros(numX+1)
PAB = np.zeros(numX+1)
thGridup=np.zeros(numX+1)
thGriddown=np.zeros(numX+1)

#a=np.loadtxt("/home/ctorney/data.txt")
#for p in np.unique(a[:,1]): print p, psw(K*np.mean(a[a[:,1]==p,2]))



l = 0.0020

lf = 2.0*(1.0-m.exp(-l*0.5))
totalP=0.0
Ed=0.0
# expected distance to NN
for dX in range(0,32):
    X = 1.0/64.0
    dist = float(dX)/numX
    
    if dist<0.5-X:
        f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
        totalP+= tup(f)
        Ed += dist*tup(f)
    elif dist<1.0-X/2.0:
        f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
        totalP+= tup(f)
        Ed += dist*tup(f)
        
print Ed/totalP
    
for fX in range(0,60):#0,numX+1):
    X = fX/float(numX)
    # find estimate of concurrent individuals

    #psolo = ((numX-fX-2.0)/float(numX))*psw((fX-nc)/float(numX-nc))
    def f(nc):
        lf = 1.0/(2.0*(1.0-m.exp(-l*0.5)))
        x2 = (fX-nc)/float(numX-nc)
        dist = nc/float(numX)
        #print x2
        xnc = x2*0.5 + lf * (1.0-x2*m.exp(-l*0.5) + (x2-1.0)* m.exp(-l*dist))
        #print xnc
        pj = (2.0/float(numX))*tup(xnc)
        psolo = ((numX-fX-2.0)/float(numX))*tup(x2)
        return nc*pj/(pj+psolo)
        #return nc - pj*fX/(pj+psolo)
    #nc1= np.max(sc.optimize.fsolve(f,fX),0)
    nc1=f(fX)
    #nc1=0
    bldist = nc1/float(numX)
    x2 = (float(fX)-nc1)/float(numX-nc1)
    #print x2#, fX, nc1
    #xnc = x2*0.5 + lf * (1.0-x2*m.exp(-l*0.5) + (x2-1.0)* m.exp(-l*dist))
    #xnc = x2
    #thGridup[fX] = (2.0/float(numX))*tup(xnc) + ((numX-fX-2.0)/float(numX))*tup(x2)
    

    ex_up=0.0
    ds = int(m.floor(0.5*(numX-nc1)))
    for d in range(1,ds):
        dist = (d+0.0)/numX
        lf = 1.0/(2.0*(1.0-m.exp(-l*0.5)))
        frac = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+bldist)))
        ex_up+=tup((1.0-frac)*x2 + frac)
    ex_up = ex_up/(ds-1.0)
    print ex_up,x2,nc1
    thGridup[fX] = (1-X)*ex_up
#        #0.5*(1.0-X)*float(av)/100.0
#    
#        if dist<0.5-X:
#            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
#            #print [X, dist, f, tup(f)]
#        elif dist<1.0-X/2.0:
#            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
#            #print [X, dist, f, tup(f)]
#        ex_up+=2.0*tup(f)
#    if (numX-fX)%2>0:
#        dist = (ds +0.0)/numX
#        if dist<0.5-X:
#            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
#            #print [X, dist, f, tup(f)]
#        elif dist<1.0-X/2.0:
#            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
#            #print [X, dist, f, tup(f)]
#        ex_up+=tup(f)
#    if X<1:
#        ex_up=ex_up/float(numX-fX)
#   
    Xs[fX] = X
#    PAB[fX] = ((f*(X)/(f*X+(1.0-X))))
#    thGridup[fX] = (1-X)*ex_up#tup(X,PAB[fX])
#    thGriddown[fX] = tdown(X,PAB[fX])
#
xGrid=np.arange(65)/64.0

plt.plot(Xs,thGridup,marker='o',label='theory up')
#pot=np.log(np.divide(thGriddown,thGridup))
#pot=np.cumsum(pot[2:64])
#plt.plot(xGrid[2:64],pot,label = str(alpha))
#plt.legend(loc=1,       ncol=1, borderaxespad=0.5,prop={'size':12})
#    
#plt.axis([0, 1, -10, 10])
#
##for p in thGridup: print p
#
##
##l = 0.0001
##
##lf = 2.0*(1.0-m.exp(-l*0.5))
##for fX in range(0,65):#0,numX+1):
##    X = fX/float(numX)
##    ex_up=0.0
##    for av in range(100):
##        dist = 0.5*(1.0-X)*float(av)/100.0
##    
##        if dist<0.5-X:
##            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
##            #print [X, dist, f]
##        elif dist<1.0-X/2.0:
##            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
##            #print [X, dist, f]
##        ex_up+=tup(f)
##    ex_up=ex_up/100.0
##    Xs[fX] = X
##    PAB[fX] = ((f*(X)/(f*X+(1.0-X))))
##    thGridup[fX] = (1-X)*ex_up#tup(X,PAB[fX])
##    thGriddown[fX] = tdown(X,PAB[fX])
##    
##plt.plot(Xs,thGridup,marker='o',label='theory up')
###plt.plot(Xs,thGriddown,marker='o',label='theory down')
###
###
###sampleP = np.load('../potential_v2/build/potential2-0.npy')
###xGrid=np.arange(65)/64.0
###yGrid = sampleP[:,0]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
####plt.plot(xGrid,yGrid,marker='o',label='sim up')
###
####yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
###yGrid = sampleP[:,1]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
####plt.plot(xGrid,yGrid,label='sim down')
###
####plt.plot(Xs,PAB)
####plt.plot(Xs,Xs)
###
###
###    
###    