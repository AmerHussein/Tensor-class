# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:11:21 2019

@author: User
"""

import scipy.integrate as integrate #efter logic commanden "as" så  kam man använda förkortningen "integrate"
import matplotlib.pyplot as plt
import numpy as np

pi = np.pi
sqrt = np.sqrt
cos = np.cos
sin = np.sin
log = np.log

g = 10
mu = 0.5
L = 2

def deriv_z(z, phi):
    u, udot = z
    return [udot, -(mu)*udot -(g/L)*u]

phi = np.linspace(0, 5.0*pi, 200000)
zinit = [pi/3.0, 0]
z = integrate.odeint(deriv_z, zinit, phi)
u, udot = z.T
fig, ax = plt.subplots()
ax.plot(phi, u)
ax.set_aspect('equal')
plt.grid(True)
plt.show()

x = np.array([1,2])
print(*x)

def maxima(a,b):
    if a>0 and b<0:
        return a
    if a<0 and b>0:
        return b
    if a>0 and b>0:
        c=a-b
        if c>0:
            return a
        if c<0:
            return b
        if a<0 and b<0:
            c=a-b
            if c>0:
                return b
            if c<0:
                return a
print(maxima(*x))

























