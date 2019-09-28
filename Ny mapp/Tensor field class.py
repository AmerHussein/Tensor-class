# -*- coding: utf-8 -*-
"""
Created on Sun Jun 30 16:46:01 2019

@author: User
"""

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import scipy.linalg as sl
from sympy import symbols, diff






















































def fastlogapprox_amer(n, x, fast_diff_mat="off"):

    a,g = [(1+x)/2], [np.sqrt(x)]

    for i in range(n):
        a_new = (a[-1] + g[-1])/2
        a.append(a_new)
        g_new = np.sqrt(g[-1]*a[-1])
        g.append(g_new)

    l = len(a)
    a = np.array(a)
    d = np.zeros((l,l))
    d[:1] = a
    for j in range(1,n+1):
        for i in range(j,n+1):
            d[j][i] = (d[j-1,i] -(d[j-1,i-1]*4**(-j)))/((1-4**(-j)))

    if fast_diff_mat == "on":
        d[-1:][0] = np.abs((x-1)/d[-1,-1] - (x-1)/a[-1])
        return d
    elif fast_diff_mat=="off":
        return (x-1)/d[-1,-1]









class Rank_1_TensorField:
    """
    Creating a class for modelling a Rank 1 tensor field, (basically a vectorfield).
    Its purpose is to capture the invraiant attribute defining a tensor. I
    Wanted to make a genereal tensor field class  but i don't know how to organise the tructure
    so i'll start with the case of rank 1 tensors. (Rank 0 tensors are normal scalars so...yeah)

    """
    def __init__(self, v=None, J=None, x=[x,y,z], B=[u,v,w], co , contra):
        """
        Here i'm trying to figure out the defining attributes of a vectorfield. I know how
        rank 1 tensors transform, with the help of the jackobian.
        """
        if not v:
            if lambda *x: self.v(*x):
                self.v = self.np.array(self.v)
                self.co = co
                self.J = lambda *x: self.J(*x)
                self.J = self.np.array(self.J)

            elif lambda *B: self.vector(*B):
                self.v = np.array(self.v)
                self.contra = contra
                self.J = lambda *y: self.J(*B)
                self.J = self.np.array(self.J)


    def phi(self):
        X = []
        X_ = [] #will contain the inverse functions of the elemts of X

        for a in x:
            a = lambda *y: a(*y)
            X.append(a)


        for b in B:
            b = lambda *x: b(*x)
            X_.append(a)

    def co_transformation(self):
        if lambda *x: self.v(*x):
            return np.matmul(self.J, self.v)
        elif lambda *B: self.vector(*B):
            return





































































































