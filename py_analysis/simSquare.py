#!/usr/bin/python

import sympy as sp
import scipy as sc
from scipy.signal import lfilter
from IPython.display import display
import numpy as np
from numpy import random
import matplotlib.pyplot as plt
from math import *
import matplotlib as mpl

from pylab import *



K = 8        
wg = 0.875
ws = 0.65


def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*math.erf(math.sqrt(wg)*(1.0-gc))
    

# Function definition is here
def tup( x ):
    return (1-x) * sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1))

def tdown( x ):
    return (x) * (1 -  sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1)))


Nx = 8
NA = Nx * Nx

NT=5000


alpha = 1




wgs = 0.01 + wg*1.0*(rand(Nx, Nx)<14.0/64.0)

a = exp(-wgs)
b = (1-a)
c = sqrt((1.0-a*a)/(2.0*wgs))



for acount in range(10,11):
    alpha = 0.0#0.1*float(acount)
    
    counts = np.zeros(NA+1)
    ups = np.zeros(NA+1)
    downs = np.zeros(NA+1)
    ccounts = np.zeros(NA+1)
    dcounts = np.zeros(NA+1)
    cups = np.zeros(NA+1)
    cdowns = np.zeros(NA+1)
    conds = np.zeros(NA+1)
    
    for av in range(100):
        print av
        G=np.ones((Nx,Nx))
        u=np.zeros((Nx,Nx))
        
        
                
                
                
        
            
        nextCount =  np.sum(u)
                        
        for t in range(1,NT):
        
            thisCount = nextCount
            G = a*G + (1.0-a) + c * randn(Nx,Nx)
            ix = randint(0,Nx)
            iy = randint(0,Nx)
            i = iy * Nx + ix
        
                    
          
            
            neigh = np.array( ( ( 1, 1 ), ( 1, 0 ), ( 1, -1 ) , ( 0, 1 ), ( 0, -1 ), ( -1, -1 ) , ( -1, 0 ), ( -1, 1 ) ))
            deltan = 0
           
            for e in range(K):
                if rand()<alpha:
                    n2 = randint(0,K)
                    
        
                    x_n = (((ix + neigh[n2,0]) % Nx) + Nx) % Nx;
                    y_n = (((iy + neigh[n2,1]) % Nx) + Nx) % Nx;
        
                    
                    if (u[x_n,y_n]>0.5):
                        deltan = deltan + 1
        
                
                else:
                    x_n = randint(0,Nx-1)
                    y_n = randint(0,Nx-1)
        
                    if (u[x_n,y_n]>0.5):
                        deltan = deltan + 1
            if  (u[ix,iy] >0.5):
                cups[thisCount] += deltan / float(K)
                ccounts[thisCount]+=1
            else:
                cdowns[thisCount] += deltan / float(K)
                dcounts[thisCount]+=1
            conds[thisCount] += deltan / float(K)
            deltan = deltan*2
            deltan = deltan - K
        
        
            pup = exp(-4.0*wgs[ix,iy]*G[ix,iy])
            pall = pup*(((1.0 - ws)/ws)**deltan)
            
            if (pall<1.0):
                u[ix,iy] = 1
            else:
                u[ix,iy] = 0
            
            counts[thisCount]+=1
        
            
            nextCount =  np.sum(u)
            if (nextCount>thisCount):
                ups[thisCount]+=1
            
            if (nextCount<thisCount):
                downs[thisCount]+=1
            
            if t % 100000000 == 0 :
                print np.sum(u)
                #plt.imshow(u, extent=[0,1,0,1], aspect='equal', vmin=0, vmax=1)
                #plt.set_cmap('hot')
                #fileName = "/home/ctorney/tmp_frames/" + '{0:05d}'.format(t) +".png"
                #plt.savefig(fileName)
    ups = np.divide(ups,counts)
    downs = np.divide(downs,counts)
    cups = np.divide(cups,ccounts)
    cdowns = np.divide(cdowns,dcounts)
    conds = np.divide(conds,counts)
    xGrid=np.arange(65)/64.0
    plt.plot(xGrid,ups,marker='o',label='sim up')
    plt.plot(xGrid,downs,marker='o',label='sim up')
    
    outdata = np.vstack((ups,downs))
    #    for r in ups: print r
    #    for r in downs: print r
    #outfile = "potential1" + '{0:02d}'.format(int(acount)) +".npy"
    #np.save(outfile, outdata)  
