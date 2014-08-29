
#!/usr/bin/python

import sympy as sp
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

a=np.loadtxt("/home/ctorney/data.txt")
for p in np.unique(a[:,1]): print p, psw(K*np.mean(a[a[:,1]==p,2]))


l = 20

lf = 2.0*(1.0-m.exp(-l*0.5))
for fX in range(32,33):#0,numX+1):
    X = fX/float(numX)
    ex_up=0.0
    ds = int(m.floor(0.5*(numX-fX)))
    for d in range(ds):
        dist = (d+0.0)/numX#0.5*(1.0-X)*float(av)/100.0
    
        if dist<0.5-X:
            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
            print [X, dist, f, tup(f)]
        elif dist<1.0-X/2.0:
            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
            print [X, dist, f, tup(f)]
        ex_up+=2.0*tup(f)
    if (numX-fX)%2>0:
        dist = (ds +0.0)/numX
        if dist<0.5-X:
            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
            print [X, dist, f, tup(f)]
        elif dist<1.0-X/2.0:
            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
            print [X, dist, f, tup(f)]
        ex_up+=tup(f)
    if X<1:
        ex_up=ex_up/float(numX-fX)
   
    Xs[fX] = X
    PAB[fX] = ((f*(X)/(f*X+(1.0-X))))
    thGridup[fX] = (1-X)*ex_up#tup(X,PAB[fX])
    thGriddown[fX] = tdown(X,PAB[fX])
    
plt.plot(Xs,thGridup,marker='o',label='theory up')
for p in thGridup: print p

#
#l = 0.0001
#
#lf = 2.0*(1.0-m.exp(-l*0.5))
#for fX in range(0,65):#0,numX+1):
#    X = fX/float(numX)
#    ex_up=0.0
#    for av in range(100):
#        dist = 0.5*(1.0-X)*float(av)/100.0
#    
#        if dist<0.5-X:
#            f = (1.0/lf)*(m.exp(-l*dist) - m.exp(-l*(dist+X)))
#            #print [X, dist, f]
#        elif dist<1.0-X/2.0:
#            f = (1.0/lf)*(m.exp(-l*dist) + m.exp(-l*(1.0-dist-X)) - 2.0*m.exp(-l*(0.5)))
#            #print [X, dist, f]
#        ex_up+=tup(f)
#    ex_up=ex_up/100.0
#    Xs[fX] = X
#    PAB[fX] = ((f*(X)/(f*X+(1.0-X))))
#    thGridup[fX] = (1-X)*ex_up#tup(X,PAB[fX])
#    thGriddown[fX] = tdown(X,PAB[fX])
#    
#plt.plot(Xs,thGridup,marker='o',label='theory up')
##plt.plot(Xs,thGriddown,marker='o',label='theory down')
##
##
##sampleP = np.load('../potential_v2/build/potential2-0.npy')
##xGrid=np.arange(65)/64.0
##yGrid = sampleP[:,0]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
###plt.plot(xGrid,yGrid,marker='o',label='sim up')
##
###yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
##yGrid = sampleP[:,1]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
###plt.plot(xGrid,yGrid,label='sim down')
##
###plt.plot(Xs,PAB)
###plt.plot(Xs,Xs)
##
##
##    
##    