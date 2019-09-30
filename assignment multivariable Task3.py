# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 21:07:03 2019

@author: User
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Sep 28 21:17:10 2019

@author: User
"""
import scipy as sc
import matplotlib.widgets  as mpw
from mpl_toolkits.mplot3d import axes3d
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from sympy import *
from sympy import symbol, diff
import time
from scipy.integrate import quad
import scipy.linalg as sl
from sympy import sqrt, sin, cos
from scipy import optimize
from scipy.linalg import qr
import matplotlib.pyplot as plot
M = np.array([[1,2,3],
              [4,5,6]])
print([*M])

def Task3(x1,x2, y1,y2, N, p):
     def z_vals(x0,y0):
          def F(x,y,z):
               return x + 2*y + z + np.exp(2*z) -1
          def g(z):
               return F(x0,y0,z)
          return optimize.fsolve(g, np.array([x0,y0]))[0]

     def taylor_coeff(F, p=p):
          def dx(f,x,y,h=1.e-8):
               return (f(x+h, y) -f(x,y))/h

          def dy(f,x,y,h=1.e-8):
               return (f(x, y+h) -f(x,y))/h

          def d2x(f,x,y,h1=1.e-8,h2=1.e-4):
               return (dx(f,x+h2,y,h1) - dx(f, x,y,h1))/h2

          def d2y(f,x,y,h1=1.e-8,h2=1.e-4):
               return (dy(f,x,y+h2,h1) - dy(f,x,y,h1))/h2

          def dxy(f,x,y,h1=1.e-8,h2=1.e-4):
               return (dx(f,x,y+h2,h1) - dx(f,x,y,h1))/h2

          def Grad(f, x,y):
               return np.array([dx(f,x,y), dy(f,x,y)])

          def Hess(f,x,y):
               H1 = np.array([d2x(f,x,y), dxy(f,x,y)])
               H2 = np.array([dxy(f,x,y), d2y(f,x,y)])
               return np.row_stack([H1, H2])

          nablaZ = Grad(F, p[0], p[1])
          hessZ = Hess(F, p[0], p[1])
          return np.column_stack([nablaZ, hessZ.T])

     def zdeg2_approx(h1,h2, p=p):
          """
          this function is the second degree polynomial approximation
          of z_vals
          """
          h    = np.array([h1,h2])
          hh   = h*h.reshape(-1,1)
          M    = taylor_coeff(z_vals)
          z    = z_vals(*p)
          zp   = np.matmul(h, M[:,0:1])[0]
          zpp  = np.trace(np.matmul(hh, M[:,1:3]))
          return z + zp + zpp

     def error(x,y):
          return np.abs(z_vals(x,y) - zdeg2_approx(x,y))/np.abs(z_vals(x,y))


     def quick_run(FUNCTION, x1,x2,y1,y2, N, p=p):
          """
          THis fucntion runs a  3d plot of the callable function FUCNTION
          """
          def limits(FUNCTION, x1=x1,x2=x2,y1=y1,y2=y2, N=N, p=p):
               x = np.linspace(x1-p[0], x2+p[0], N)
               y = np.linspace(y1-p[1], y2+p[1], N)
               X,Y  = np.meshgrid(x,y)
               Z = np.array([[FUNCTION(X[i,j]-p[0], Y[i,j]-p[1]) for i in range(X.shape[0])] for j in range(Y.shape[0])])
               return X, Y, Z

          fig = plt.figure()
          ax = fig.gca(projection="3d")
          for func in FUNCTION:
               ax.plot_surface(limits(func)[0], limits(func)[1], limits(func)[2], alpha=0.8)

          #set limits
          ax.set_xlim3d(x1,x2)
          ax.set_ylim3d(y1,y2)

          #set labels
          ax.set_xlabel("X axis")
          ax.set_ylabel("Y axis")
          ax.set_zlabel("Z axis")
          plt.show()
     return quick_run([zdeg2_approx, z_vals, error],-1,1,-1,1,100)

print(Task3(-1,1,-1,1,100,np.array([0,0.])))

















