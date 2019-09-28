# -*- coding: utf-8 -*-
"""
Created on Mon May 20 09:26:22 2019

@author: User
"""

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np


#x1 = np.linspace(1, 100, 10000)
#x2 = np.linspace(-50, -1, 10000)

#figure()

#subplot(211)
#plot(x1, (1+x1**2)/x1/(1+x1)**2, "o")
#plot(x2, (1+x2**2)/x2/(1+x2)**2, ".")
#subplot(212)



#legend(loc = "center" , fontsize = "small")

import numpy as np
import matplotlib.pyplot as plt

max_iter = 10000
step = 0.001

x = 4*np.random.random((10,2))

plt.xlim = [-9000, 10000]

def f(x, y):
    return np.array((x**3-3*x*y**2-1, 3*x**2*y-y**3))

def J(x, y):
    return np.array(((3*x**2-3*y**2, -6*x*y),
            (6*x*y, 3*x**2-3*y**2)))

class fractal2D:
    def __init__(self, F, J):
        self.F = F
        self.J = J
        self.zeros = []

    def newton(self, x0, tol=1.e-4, step=0.001, max_iter=10000):
        x = np.array(x0, dtype=float)
        F_val = self.F(*x)
        F_norm = np.linalg.norm(F_val)
        k = 0
        while abs(F_norm) > tol and k < max_iter:
            try:
                delta = np.linalg.solve(self.J(*x), -F_val)
            except np.linalg.LinAlgError:
                y = x + step
                while not np.linalg.det(self.J(*y)):
                    y += step
                delta = np.linalg.solve(self.J(*y), -F_val)
            x += delta
            F_val = self.F(*x)
            F_norm = np.linalg.norm(F_val)
            k += 1
        if abs(F_norm) > tol:
            return np.array(np.inf)
        return x

    def zero_index(self, x0, tol=1.e-4):
        x = self.newton(x0)
        if (x == np.inf).any():
            return np.inf
        if not len(self.zeros):
            self.zeros.append(x)
            return 0
        for k,i in enumerate(self.zeros):
            if np.linalg.norm(i-x) < tol:
                return k
        self.zeros.append(x)
        return len(self.zeros) - 1

    def plot(self, N, a, b, c, d):
        x = np.linspace(a, b, N)
        y = np.linspace(c, d, N)
        X, Y = np.meshgrid(x, y)
        A = np.array([np.array([self.zero_index(np.array([X[i,j],Y[i,j]])) for i in range(N)]) for j in range(N)])
        fig, ax = plt.subplots()
        ax.pcolor(X, Y, A)
frac = fractal2D(f, J)
frac.plot(1000, *(-1,1), *(-1,1))

plt.show()


































































































































































