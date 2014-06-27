
#!/usr/bin/python

import sympy as sp
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

K = 8        
wg = 0.2875
ws = 0.65
alpha = 1



for acount in range(0,11):
    alpha = 0.1*float(acount)
    

    sampleP = np.load('potential' + '{0:02d}'.format(int(10.0*alpha)) + '.npy')
    xGrid=np.arange(65)/64.0
    
    
    pot=np.log(np.divide(sampleP[1,:],sampleP[0,:]))
    pot=np.cumsum(pot[2:64])
    plt.plot(xGrid[2:64],pot)
    potmatch = pot[0]
    
    plt.axis([0, 1, -10, 10])
    

def tup( x ):
    alpha1 = alpha*1.0/(1.0+exp(-(x-xthresh)/thresh))
    if x > 0.0625:
        alpha1=alpha
    else:
        alpha1 = 0
    alpha1 = alpha*(1.0-exp(-(x)/0.1))
    xx = (1.0-alpha1)*(x)
    return (1-x) * sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))

def tdown( x ):
    alpha1 = alpha*1.0/(1.0+exp(-(x-xthresh)/thresh))
    if x > 0.0625:
        alpha1=alpha
    else:
        alpha1 = 0
    alpha1 = alpha*(1.0-exp(-(x)/0.1))        
    xx = 1.0-(1.0-alpha1)*(1.0-x)
    return (x) * (1 -  sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1)))


figure()
for acount in range(0,6):
    alpha = 0.15*float(acount)
    
    for x in range(0,65):
        thGridup[x] = tup(x/64.0)
        thGriddown[x] = tdown(x/64.0)

    xGrid=np.arange(65)/64.0
    
    
    pot=np.log(np.divide(sampleP[1,:],sampleP[0,:]))
    pot=np.log(np.divide(thGriddown,thGridup))
    pot=np.cumsum(pot[2:64])
    plt.plot(xGrid[2:64],pot)
    potmatch = pot[0]
    
    plt.axis([0, 1, -10, 10])

