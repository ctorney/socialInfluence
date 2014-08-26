
#!/usr/bin/python

import sympy as sp
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

K = 8        
wgb = 0.2
ws = 0.65
alpha = 0.0
NA = 64
def psw( j, wggg):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wggg)
    return 0.5 + 0.5*math.erf(math.sqrt(wggg)*(1.0-gc))
    

# Function definition is here
def tup( x ):
    alpha1 = alpha*(1.0-exp(-x/0.1))
    xx = (1.0-alpha1)*(x)
    return 0.5*(1-x) * sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j,wgb) for j in xrange(0,K+1)) +0.5*(1-x) * sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j,wg) for j in xrange(0,K+1))

def tdown( x ):
    alpha1 = alpha*(1.0-exp(-x/0.1))        
    xx = 1.0-(1.0-alpha1)*(1.0-x)
    return 0.5*(x) * (1 -  sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j,wgb) for j in xrange(0,K+1))) + 0.5*(x) * (1 -  sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j,wg) for j in xrange(0,K+1)))



numA = 6
cost = 000000000
alphas = np.zeros(numA)
atimes = np.zeros(numA)
for acount in range(0,numA):
    alpha = float(acount)/5.0
    numWG = 20
    atimes = np.zeros(numWG)
    wgs = np.zeros(numWG)
    for wgi in range(numWG):

        wg = wgb + wgi*0.5/float(numWG)
        N2 = int(NA*0.5+1)
    
        Q = np.zeros((N2,N2))

        Q[0,0]=1.0-tup(0)-tdown(0)
        Q[0,1]=tup(0)

        Q[N2-1,N2-2]=tdown(N2/float(NA))
        Q[N2-1,N2-1]=1.0-tup(N2/float(NA))-tdown(N2/float(NA))
        for x in range(1,N2-1):
            Q[x,x-1] = tdown(x/float(NA))
            Q[x,x] = 1.0 - tup(x/float(NA)) - tdown(x/float(NA))
            Q[x,x+1] = tup(x/float(NA))
    
        bb = np.matrix(np.linalg.inv(np.identity(N2) - Q))

    
        times = bb*np.matrix(np.ones((N2,1)))


        atimes[wgi] = times[0] + cost*wgi*0.5/float(numWG)
        #print times[0]
        wgs[wgi] = wg - wgb
    f = plot(wgs,atimes,label = str(alpha))
    plt.yscale('log')
    legend(loc=1,       ncol=1, borderaxespad=0.5,prop={'size':12})

axis([0, 0.31, 100, 10000000000000])
plt.savefig('fig1.eps', format='eps', dpi=100)