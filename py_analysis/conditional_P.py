
#!/usr/bin/python

import sympy as sp
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

K = 8        
wg = 0.5
ws = 0.65
alpha = 0.0
NA = 64.0
N = 1024
def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*math.erf(math.sqrt(wg)*(1.0-gc))
    

# Function definition is here
def tup( x ):
    alpha1 = alpha*(1.0-exp(-x/0.1))
    xx = (1.0-alpha1)*(x)
    return (1-x) * sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))

def tdown( x ):
    alpha1 = alpha*(1.0-exp(-x/0.1))        
    xx = 1.0-(1.0-alpha1)*(1.0-x)
    return (x) * (1 -  sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1)))

def pswB1( x ):
    return  (sum(sp.binomial(K-1,j) * x**j * (1-x)**(K-1-j) * psw(j+1) for j in xrange(0,K)))
def pswNB1( x ):
    return  (sum(sp.binomial(K-1,j) * x**j * (1-x)**(K-1-j) * psw(j) for j in xrange(0,K)))
    
def pswB( x ):
    xx=alpha*(1.0+(x-1.0/NA)*(float(K)-1.0))/float(K) + (1-alpha)*(x-1.0/NA)
    return sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))
def pswNB( x ):
    xx=alpha*(x*(float(K)-1.0))/float(K) + (1-alpha)*x
    return sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))
def pswBNB( x ):
    return sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1))

numX = 100
Xs = np.zeros(numX+1)
PAB = np.zeros(numX+1)
for fX in range(0,numX+1):
    X = fX/float(numX)
    alpha = 1
    pInt = alpha
    f = (pInt*pswB(X) + (1.0-pInt)*pswBNB(X))/(pInt*pswNB(X) + (1.0-pInt)*pswBNB(X))
    Xs[fX] = X
    PAB[fX] = ((f*(X)/(f*X+(1.0-X)))-X)/((1-X))
    
plt.plot(Xs,PAB)
#plt.plot(Xs,Xs)


    
    