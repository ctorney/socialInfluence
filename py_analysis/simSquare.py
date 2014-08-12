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
wg = 0.5
ws = 0.75


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


NT=4000  


alpha = 1

counts = np.zeros(NA+1)
ups = np.zeros(NA+1)
downs = np.zeros(NA+1)




a = exp(-wg)
b = (1-a)
c = sqrt((1.0-a*a)/(2.0*wg))




for av in range(50):
        
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
     
        deltan = deltan*2
        deltan = deltan - K
    
    
        pup = exp(-4.0*wg*G[ix,iy])
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

xGrid=np.arange(65)/64.0
plt.plot(xGrid,ups,label='ups')


plt.plot(xGrid,downs,label='downs')
plt.axis([0, 0.5, 0, 0.500])
legend()
figure()
pot=np.log(np.divide(downs,ups))
pot=np.cumsum(pot[2:64])


plt.plot(xGrid[2:64],pot)

thGridup=np.zeros(65)
thGriddown=np.zeros(65)

for x in range(0,65):
    thGridup[x] = tup(x/64.0)
    thGriddown[x] = tdown(x/64.0)
#plt.plot(xGrid,thGridup)
#plt.plot(xGrid,thGriddown)

pot=np.log(np.divide(thGriddown,thGridup))
pot=np.cumsum(pot[2:64])
plt.plot(xGrid[2:64],pot)