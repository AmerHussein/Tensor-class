# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 12:23:47 2019

@author: User
"""





from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff


class Interval:
    def __init__(self, l_b, u_b):
        self.l_b = l_b
        self.u_b = u_b
      #  if l_b or u_b is None:
      #      self.l_b = float(l_b)
       #     self.u_b = float(u_b)

    def __max_min__(a,b, shift):
        if a>0 and b<0 and shift==1:
            return a
        elif a>0 and b<0 and shift==-1:
            return b
        if a<0 and b>0 and shift==1:
            return b
        elif a<0 and b>0 and shift==-1:
            return a
        if a>0 and b>0 and shift==1:
            c=a-b
            if c>0:
                return a
            if c<0:
                return b
        elif a>0 and b>0 and shift==-1:
            c=a-b
            if c>0:
                return b
            if c<0:
                return a
        if a<0 and b<0 and shift==1:
            c=-(a-b)
            if c>0:
                return b
            if c<0:
                return a
        elif a<0 and b<0 and shift==-1:
            c=a-b
            if c>0:
                return a
            if c<0:
                return b

    def __repr__ (self):
        M = [self.l_b, self.u_b]
        return "[{} , {}]".format(M[0], M[1])


    def __add__(self, other):
        M = [self.l_b, self.u_b]
        N = [other.l_b, other.u_b]
        return Interval(M[0]+N[0], M[1]+N[1])

    def __sub__(self, other):
        M = [self.l_b, self.u_b]
        N = [other.l_b, other.u_b]
        return Interval(M[0]-N[1],M[1]-N[0])

    def __mult__(self, other):
        M = [self.l_b, self.u_b]
        N = [other.l_b, other.u_b]
        Q = np.zeros((2,2))
        for i in range(2):
            for j in range(2):
                Q[i][j] = M[i]*N[j]

        r = np.array([ max_min(*Q[1],1), max_min(*Q[0],-1)])
        return Interval(*r)

    def __truediv__(self, other, tol=10**(-4)):
        M = np.array([self.l_b, self.u_b])
        N = np.array([other.l_b, other.u_b])
        m1 = np.array([{},{}])
        m2 = np.array([{},{}])
        m12 = np.array([{},{}])
        R = np.zeros((2,2))

       # for k in range(1):
        #    s[k] = abs(t[k] - N[k])
        #if s[0] < tol and s[1] < tol:
        for i in range(2):
            for j in range(2):
                R[i][j] = M[i]/N[j]

        for k in range(2):
            m1[k] = max_min(*R[k],-1)
            m2[k] = max_min(*R[k], 1)
        m12 = [max_min(*m2, 1), max_min(*m1,-1)]

      #  r =  np.array([max_min(*[max_min(*R[0],-1), max_min(*R[1],-1)],-1),
       #                 max_min(*[max_min(*R[0],1), max_min(*R[1],1)],1)])

       # else:
        #    if s[0] < tol or s[1] < tol:
         #       return np.array((np.inf,np.inf))
        return Interval(*m12), R, m12

    def __contains__(other, x):
        N = np.array([other.l_b, other.u_b])
        if x == N[0] or x ==N[1]:
            if x < N[1]:
                return "number is the lower bound"
            if x > N[0]:
                return "number is the upper bound"

        if x !=N[0] and x !=N[1]:
            delta1 = x-N[0]
            if delta1 <= 0:
                return"Number is too low"
            elif delta1 >=0:
                delta2 = N[1]-x
                if delta2 <=0:
                    return "number is to high"
                if delta2 >=0:
                    return "Interval Contains number"
                    return x

    def __Sensordiv__(self, other):
        tol = 1/10**4
        M = np.array([self.l_b == float, self.u_b==float])
        N = np.array([other.l_b == float, other.u_b== float])
        m1 = np.array([{},{}])
        m2 = np.array([{},{}])
        m12 = np.array([{},{}])
        R = np.zeros((2,2))

        #if abs(N[0]) < tol or abs(N[1]) < tol:
        #    raise ZeroDivisionError("can't divide by Zero")

        for i in range(2):
            for j in range(2):
                R[i][j] = M[i]/N[j]
        for k in range(2):
            m1[k] = max_min(*R[k],-1)
            m2[k] = max_min(*R[k], 1)
        m12 = [max_min(*m2, 1), max_min(*m1,-1)]

        if abs(m12[0])> 1/tol and abs(m12[1])<1/tol:
            raise OverflowError("lower bound diverges")
            return Interval(np.inf, m12[1])

        elif abs(m12[0])<1/tol and abs(m12[1])>1/tol:
            raise OverflowError("upper bound diverges")
            return Interval(m12[0], np.inf)
        elif abs(m12[0])<1/tol and abs(m12[1])<1/tol:
            return Interval(*m12), R, m12





I1 = Interval(1, 4)
I2 = Interval(-2, -1)
I3 = Interval.__add__(I1,I2)
I4 = Interval(1/2, 4)
I5 = Interval(3, 5)
print(I1, I2, I3)
print(Interval.__add__(I1, I2), Interval.__sub__(I1, I2), Interval.__mult__(I1, I2), Interval.__truediv__(I1,I2))
x = 4
print(Interval.__contains__(I1, x))
print(Interval.__Sensordiv__(I5, I4))
































