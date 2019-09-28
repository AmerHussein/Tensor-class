# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 22:32:16 2019

@author: User
"""



from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np




def fastlogapprox(x, n):
    a = np.zeros(n)
    g = np.zeros(n)
    d = np.zeros((n,n))
    a[0] = 0.5*(1+x[1])
    g[0] = np.sqrt(x[1])
    for i in range(1, n):
        a[i] = 0.5*(a[i-1]+g[i-1])
        g[i] = np.sqrt(a[i]*g[i-1])
        d[0:,i] = a
    for i in range(1,n):
        d[i,i] = (a[i] - 4**(-i)*a[i-1])/(1-4**(-i))
    return (x-1)/d[-1,-1]




x_vals = np.linspace(0, 20, 1000)
for k in range(1, 5):
    y_vals = np.abs(fastlogapprox(x_vals, k) - np.log(x_vals))
    plt.plot(x_vals, y_vals)
    plt.title('Error behavior of the accelerated Carlsson Method for the log')
    plt.xlabel('x')
    plt.ylabel('error')
    plt.legend(loc='upper left')
    plt.show()










































