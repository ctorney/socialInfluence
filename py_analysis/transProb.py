
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
alpha = 0.00000450
thresh = 0.125
xthresh = 0.06125
NA = 64
def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*math.erf(math.sqrt(wg)*(1.0-gc))
    

# Function definition is here
def tup( x ):
    alpha1 = alpha*1.0/(1.0+exp(-(x-xthresh)/thresh))
    if x > 0.0625:
        alpha1=alpha
    else:
        alpha1 = 0
    alpha1 = alpha*(1.0-exp(-x/0.1))
    xx = (1.0-alpha1)*(x)
    return (1-x) * sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))

def tdown( x ):
    alpha1 = alpha*1.0/(1.0+exp(-(x-xthresh)/thresh))
    if x > 0.0625:
        alpha1=alpha
    else:
        alpha1 = 0
    alpha1 = alpha*(1.0-exp(-x/0.1))        
    xx = 1.0-(1.0-alpha1)*(1.0-x)
    return (x) * (1 -  sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1)))


sampleP = np.load('potential10.npy')
xGrid=np.arange(65)/64.0
yGrid = sampleP[0,:]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
plt.plot(xGrid,yGrid,label='sim up')

#yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
yGrid = sampleP[1,:]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
plt.plot(xGrid,yGrid,label='sim down')
plt.axis([0, 1, 0, 0.4])

sampleP = np.load('potential00.npy')
xGrid=np.arange(65)/64.0
yGrid = sampleP[0,:]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
plt.plot(xGrid,yGrid,label='sim up')

#yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
yGrid = sampleP[1,:]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
plt.plot(xGrid,yGrid,label='sim down')
plt.axis([0, 1, 0, 0.4])


thGridup=np.zeros(65)
thGriddown=np.zeros(65)

for x in range(0,65):
    thGridup[x] = tup(x/64.0)
    thGriddown[x] = tdown(x/64.0)
plt.plot(xGrid,thGridup,label='theory up')
plt.plot(xGrid,thGriddown,label='theory down')
legend()
plt.figure()

sampleP = np.load('potential00.npy')
pot=np.log(np.divide(sampleP[1,:],sampleP[0,:]))
pot=np.cumsum(pot[2:64])
plt.plot(xGrid[2:64],pot)

potmatch = pot[0]
#pot=np.log(np.divide(sampleP[1,:],sampleP[0,:]))
pot=np.log(np.divide(thGriddown,thGridup))
pot=np.cumsum(pot[2:64])
potdiff = potmatch-pot[0]
pot = pot + potdiff
plt.plot(xGrid[2:64],pot)

plt.axis([0, 1, -4, 10])
