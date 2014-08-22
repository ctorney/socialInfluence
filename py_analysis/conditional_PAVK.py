
#!/usr/bin/python

import sympy as sp
from IPython.display import display
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

K = 80       
wg = 0.25
ws = 0.51547
alpha = 0.0
NA = 64.0
N = 1024
def psw( j ):
    gc = np.log(ws/(1-ws))*(K-2*j)/(4*wg)
    return 0.5 + 0.5*math.erf(math.sqrt(wg)*(1.0-gc))
    

# Function definition is here
def tup( x,xx ):
    #alpha1 = alpha*(1.0-exp(-x/0.1))
    #xx = (1.0-alpha1)*(x)
    xxx=x
    if x<1:
        xxx = (1.0-xx)*(x)/(1.0-x)
    
    #xxx = (1.0-alpha1)*(x)
    
    return (1-x) * psw(K*xxx)

def tdown( x,xx ):
    #alpha1 = alpha*(1.0-exp(-x/0.1))        
    #xx = 1.0-(1.0-alpha1)*(1.0-x)
    return (x) * (1 -  psw(K*xx))# sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1)))

def pswB1( x ):
    return  (sum(sp.binomial(K-1,j) * x**j * (1-x)**(K-1-j) * psw(j+1) for j in xrange(0,K)))
def pswNB1( x ):
    return  (sum(sp.binomial(K-1,j) * x**j * (1-x)**(K-1-j) * psw(j) for j in xrange(0,K)))
    
def pswB( x ):
    xx=alpha*(1.0+(x)*(float(K)-1.0))/float(K) + (1-alpha)*(x)
    return sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))
def pswNB( x ):
    xx=alpha*(x*(float(K)-1.0))/float(K) + (1-alpha)*x
    return sum(sp.binomial(K,j) * xx**j * (1-xx)**(K-j) * psw(j) for j in xrange(0,K+1))
def pswBNB( x ):
    return sum(sp.binomial(K,j) * x**j * (1-x)**(K-j) * psw(j) for j in xrange(0,K+1))

numX = 64
Xs = np.zeros(numX+1)
PAB = np.zeros(numX+1)
thGridup=np.zeros(numX+1)
thGriddown=np.zeros(numX+1)


a=[0,0.00151943,0.031433,0.0639558,0.0968939,0.124662,0.158386,0.181953,0.215152,0.245279,0.265507,0.291785,0.303819,0.345529,0.347388,0.36632,0.387305,0.378803,0.428597,0.396579,0.448437,0.464057,0.49104,0.491688,0.505208,0.507632,0.530716,0.551768,0.559821,0.562003,0.567448,0.596237,0.642578,0.651515,0.632353,0.659286,0.68287,0.665848,0.733553,0.674679,0.707102,0.73374,0.720588,0.739535,0.754261,0.770833,0.783062,0.764184,0.796131,0.798469,0.825,0.834069,0.839103,0.862264,0.870581,0.879545,0.898252,0.899342,0.916188,0.917373,0.936905,0.953454,0.970033,0.984628,1]


for fX in range(0,numX+1):
    X = fX/float(numX)
    alpha = 0
    pInt = alpha
    f = (pInt*pswB(X) + (1.0-pInt)*pswBNB(X))/(pInt*pswNB(X) + (1.0-pInt)*pswBNB(X))
    Xs[fX] = X
    PAB[fX] = ((f*(X)/(f*X+(1.0-X))))
    thGridup[fX] = tup(X,PAB[fX])
    thGriddown[fX] = tdown(X,PAB[fX])
    
plt.plot(Xs,thGridup,marker='o',label='theory up')
plt.plot(Xs,thGriddown,marker='o',label='theory down')


sampleP = np.load('../potential_v2/build/potential2-0.npy')
xGrid=np.arange(65)/64.0
yGrid = sampleP[:,0]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
#plt.plot(xGrid,yGrid,marker='o',label='sim up')

#yGrid = np.concatenate((sampleP[0,1:64:2],sampleP[1,0:64:2]))
yGrid = sampleP[:,1]#np.concatenate((sampleP[0,0:64:2],sampleP[1,1:64:2]))
#plt.plot(xGrid,yGrid,label='sim down')

#plt.plot(Xs,PAB)
#plt.plot(Xs,Xs)


    
    