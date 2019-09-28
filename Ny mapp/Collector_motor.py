# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 14:24:17 2019

@author: User
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 22:42:29 2019

@author: User
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:50:35 2019

@author: User
"""

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




















"""
def Test_F3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([ x**8 -28*x**6*y**2 + 70*x**4*y**4 +15*x**4 -28*x**2*y**6 -90*x**2*y**6 + y**8 + 15*y**4 -16,  8*x**7*y -56*x**5*y**3 + 56*x**3*y - 8*x*y**7 -60*x*y**3])

def Test_J3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[8*x**7 - 168*x**5*y**2 + 280*x**3*y**4 + 60*x**3 - 236*x*y**6, -56*x**6*y + 280*x**4*y**3 - 708*x**2*y**5 + 8*y**7 + 60*y**3 ],
                     [56*x**6*y - 280*x**4*y**3 + 168*x**2*y - 8*y**7 - 60*y**3, 8*x**7 - 168*x**5*y**2 + 56*x**3 - 56*x*y**6 - 180*x*y**2]])


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])


def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])



print(Test_F3(np.array([ 1.00000010e+00, -7.79084175e-07])))


"""
l = [2,3]
print(l.index(2))


def Newt_zero(F, J, a,b,c,d,N,v=None, Array=[], Ind=[], option=None, margin=1.e-4, max_iter=1000, PLOT="both", zero=[]):
    def Newt_motor(F, J, v0, max_iter=max_iter, margin=margin, step=0.0001, h=1.e-9, finite=False, once=False):
        v0 = np.array(v0)
        xvals = [v0]
        V = xvals[-1]
        V_nrm = np.linalg.norm(xvals[-1])
        F_nrm = np.linalg.norm(F(xvals[-1]))
        x_max = 1/margin
        F_max = 1/margin
        k = 0
        def solving_motor(v):
            def finite_difference(v, h):
                J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                              [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                return J
            if finite == False:
                if once==False:
                    delta = np.linalg.solve(J(v), -F(v))
                    new = delta + v
                    return new
                elif once==True:
                    J_once = J(v0)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

            elif finite == True:
                if once ==False:
                    delta = np.linalg.solve(finite_difference(v, h), -F(v))
                    new = delta + v
                    return new
                elif once == True:
                    J_once = finite_difference(v0, h)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

        while k < max_iter and np.linalg.norm(xvals[-1]) < x_max and np.linalg.norm(F(xvals[-1])) > margin:
            try:
                new = solving_motor(xvals[-1])
                xvals.append(new)

            except sl.LinAlgError:
                y = xvals[-1] + step
                xvals.append(y)
            k += 1

        if np.linalg.norm(F(xvals[-1])) > margin or np.linalg.norm(xvals[-1]) > x_max:
            return (np.array((np.inf, np.inf)), k)

        else:
            return xvals[-1], k

    if not a and b and c and d and N and v.ndim==1:
        new = Newt_motor(F, J, v)
        v_a = new[0]
        v_i = new[1]
        Array.append(v_a)
        Ind.append(v_i)

        if option=="Array":
            return Array
        elif  option=="Ind":
            return Ind

    elif not v:
        def zero_collect_motor(F, J, V):
            v_a, v_i = Newt_motor(F, J, v)
            if (v_a==np.inf).any():
                return (np.inf, v_i)
            if not len(zero):
                zero.append(v_a)
                return (0,v_i)
            for ind, s in enumerate(zero):
                if np.linalg.norm(s-v_a) < margin:
                    return (ind, v_i)
            zero.append(v_a)
            return (len(zero)-1, v_i)


        def tensor_motor(a, b, c, d, N):
            x = np.linspace(a, b, N)
            y = np.linspace(c, d, N)
            L = []
            X, Y = np.meshgrid(x, y)

            for i in range(N):
                L.append([])
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    L[i].append((Y[i, j], X[i, j]))
            L = np.array(L)
            return L

        v = tensor_motor(a,b,c,d,N)
        Jesus = np.array([zero_collect_motor(F,J, v[i,j])[0] for i in range(v.shape[0])] for j in range(v.shape[1]))

        if PLOT == True:
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, Jesus, cmap='plasma')

        elif PLOT == False:
            return Jesus

        elif PLOT =="both":
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, Jesus, cmap='plasma')

            return Jesus







def Test_F3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([ x**8 -28*x**6*y**2 + 70*x**4*y**4 +15*x**4 -28*x**2*y**6 -90*x**2*y**6 + y**8 + 15*y**4 -16,  8*x**7*y -56*x**5*y**3 + 56*x**3*y - 8*x*y**7 -60*x*y**3])

def Test_J3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[8*x**7 - 168*x**5*y**2 + 280*x**3*y**4 + 60*x**3 - 236*x*y**6, -56*x**6*y + 280*x**4*y**3 - 708*x**2*y**5 + 8*y**7 + 60*y**3 ],
                     [56*x**6*y - 280*x**4*y**3 + 168*x**2*y - 8*y**7 - 60*y**3, 8*x**7 - 168*x**5*y**2 + 56*x**3 - 56*x*y**6 - 180*x*y**2]])


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])


def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])
print(Newt_zero(Test_F3, Test_J3, -1, 1, -1, 1, 10))




























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





def collect_insert_motor(F, J, a,b,c,d,N, v=None,margin=1.e-4, max_iter=1000, Array=[]):
    x_max = 1/margin
    def Newt_motor(F, J, v0, max_iter=max_iter, margin=margin, step=0.0001, h=1.e-9, finite=False, once=False):
        v0 = np.array(v0)
        xvals = [v0]
        V = xvals[-1]
        V_nrm = np.linalg.norm(xvals[-1])
        F_nrm = np.linalg.norm(F(xvals[-1]))
        x_max = 1/margin
        F_max = 1/margin
        k = 0
        def solving_motor(v):
            def finite_difference(v, h):
                J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                              [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                return J
            if finite == False:
                if once==False:
                    delta = np.linalg.solve(J(v), -F(v))
                    new = delta + v
                    return new
                elif once==True:
                    J_once = J(v0)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

            elif finite == True:
                if once ==False:
                    delta = np.linalg.solve(finite_difference(v, h), -F(v))
                    new = delta + v
                    return new
                elif once == True:
                    J_once = finite_difference(v0, h)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

        while k < max_iter and np.linalg.norm(xvals[-1]) < x_max and np.linalg.norm(F(xvals[-1])) > margin:
            try:
                new = solving_motor(xvals[-1])
                xvals.append(new)

            except sl.LinAlgError:
                y = xvals[-1] + step
                xvals.append(y)
            k += 1

        if np.linalg.norm(F(xvals[-1])) > margin or np.linalg.norm(xvals[-1]) > x_max:
            return (np.array((np.inf, np.inf)), k)

        else:
            return xvals[-1], k

    x_max = 1/margin
    if not a and b and c and d and N and v.ndim == 1:
        new = Newt_motor(F,J,v)
        v_array = new[0]
        v_ind   = new[1]
        if np.linalg.norm(F(v_array)) < margin:
            return v_array
        elif np.linalg.norm(F(v_array)) > margin or np.linalg.norm(v_array) > x_max:
            raise Exception(f"{v} does not converge to any zero")

    else:
#----------------Motors-----------------------------------------------
        def tensor_motor(a, b, c, d, N):
            x = np.linspace(a, b, N)
            y = np.linspace(c, d, N)
            L = []
            X, Y = np.meshgrid(x, y)

            for i in range(N):
                L.append([])
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    L[i].append((Y[i, j], X[i, j]))
            L = np.array(L)
            return L

        def for_motor(A, start1=0, start2=0):
            matrix = np.array((A.shape[0], A.shape[1]))
            for i in range(start1,A.shape[0]):
                for j in range(start2, A.shape[1]):
                    new = Newt_motor(F, J, A[i,j])
                    n_abs = np.linalg.norm(new[0])
                    for s in Array:
                        if np.abs(s - n_abs) < margin:
                            matrix[i,j] = x_max
                        else:
                            matrix[i,j] = n_abs
            return matrix
        def gauge_norm_motor(vector):
                first = Newt_motor(F, J, vector)
                f_a = first[0]
                f_i = first[1]
                if np.linalg.norm(F(f_a)) < margin:
                    return "True", f_a
                else:
                    return "False"
#-------------------------------------------------------------------
        e = tensor_motor(a,b,c,d,N)
        var1 = gauge_norm_motor(e[0,0])
        if var1[0]=="True":
            f_abs = np.linalg.norm(var1[1])
            Array.append(f_abs)
        else:
            for i in range(1, e.shape[0]):
                for j in range(1, e.shape[1]):
                    var2 = gauge_norm_motor(e[i,j])
                    if var2[0] =="True":
                        f_abs = np.linalg.norm(var1[1])
                        Array.append(f_abs)
            if len(Array) == 0:
                raise Exception("the meshgrid contains no zero at all,the list; Array, is empty")
            else:
                M = for_motor(e)
                return M




def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])


def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])




print(collect_insert_motor(Test_F1, Test_J1, 10,10,10,10,2))
