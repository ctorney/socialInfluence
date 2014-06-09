
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

def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*math.erf(math.sqrt(wg)*(1.0-gc))
    

# Function definition is here
def tup( x ):
    return (1-x) * sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1))

def tdown( x ):
    return (x) * (1 -  sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1)))


sampleP = np.load('potential5-50.npy')
xGrid=np.arange(64)/64.0
yGrid = np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
plt.plot(xGrid,yGrid)

yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
plt.plot(xGrid,yGrid)
plt.axis([0, 0.2, 0, 0.100])

thGridup=np.zeros(64)
thGriddown=np.zeros(64)

for x in range(0,64):
    thGridup[x] = tup(x/64.0)
    thGriddown[x] = tdown(x/64.0)
plt.plot(xGrid,thGridup)
plt.plot(xGrid,thGriddown)



