# -*- coding: utf-8 -*-
"""
Created on Sun May 26 19:59:41 2019

@author: User
"""



from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff


class CTM:
    def __init__(self, X, B, C):
         if not isinstance(B, int) and not isinstance(B, float):
             raise TypeError("the first variable should be a number")
         if not isinstance(C, int) and not isinstance(C, float):
            raise TypeError("the second variable should be a number")
         if not hasattr(X, "__len__"):
             raise TypeError("X should be a vector")

         self.X = X, self.B = B, self.C = C

    def __f__(self):
        X = zeros(1)
        def x(B, C):
            self.X[0] = x(B, C)
            return self.X[0]
        def y(B,C):
            self.x[1] = y(B,C)
            return self.X[1]
        return self.X

    def __J__(self):
        B, C = symbols("B C", real = True)
        J = np.array([diff(CTM(self.X)[0], B), diff(CTM(self.X)[0], C)]
                     [diff(CTM(self.X)[1], B), diff(CTM(self.X)[1], C)])
        return J


class Map:
    def __init__(self, f1, f2):
        self.f1 = f1
        self.f2 = f2

    def __repr__ (self):
        return "[{}, {}]".format(self.f1, self.f2)


    def __Jac__(self, B,C):
        B, C = symbols("B C", real = True)
        def Df1(B,C):
            return [diff(self.f1, B), diff(self.f1, C)]
        def Df2(B,C):
            return [diff(self.f2, B), diff(self.f2, C)]
        return np.array([Df1(self.f1), Df2(self.f2)])

class FJ:
    def __init__(self, F, J):
        if not hasattr(F, "__len__"):
            raise TypeError("F should be a vector")
        if not hasattr(J, "__shape__"):
            raise TypeError("J should be a matrix")
        self.F = np.array(F)
        self.J = np.array(J)

        def __F_to_J__(self):
            self.J  = zeros((2,2))
            self.J[0] = self.f[]






































































































































































































































