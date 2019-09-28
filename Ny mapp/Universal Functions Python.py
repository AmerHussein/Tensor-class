# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 19:20:41 2019

@author: User
"""

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff
import scipy.linalg as sl













def heaviside(x):
    if x >= 0:
        return 1.
    else:
        return 0.
vheaviside = np.vectorize(heaviside)
Q = vheaviside(np.array([-1, 2]))
print(Q)




def find_index(A,B, error_margin=0.5):
    def vfind_value(find_value):
        def find_value(x):
            if x<error_margin:
                return x
        return np.vectorize(find_value)
    q = vfind_value(A)
    return q
#    q = vfind_value(A)
#    Q = list(q)
#    _min_ = min(Q)
#    ind = A.index(_min_)
#    zero = B[ind]
#    return zero

A = np.array([5,4,3,2,1,0,1,2,3,4,5,6])
B = np.array([1,2,3,4,5,6,7,8,9,10,11,12])
print(find_index(A,B))

def zero(f, a,b, div= 1000000, error_margin = 10**(-6)):
    vf = np.vectorize(f)
    x = np.linspace(a,b,div)
    f_vals = vf(x)

    def find_lowest_value(x):
        if x < error_margin:
            return np.abs(x)
        else:
            return np.inf

    vfind_lowest_value = np.vectorize(find_lowest_value)
    Q = vfind_lowest_value(f_vals)
    q = list(Q)
    _min_ = min(q)
    ind = q.index(_min_)
    return x[ind]


def p(x):
    return 3*x**2 -5

print(zero(p,0,4))








def find_value(x, error_margin=2):
    if x < error_margin:
        return x
    else:
        return np.inf

























def vv(A, B):
    if type(A) != type(B) != np.array():
        A = np.array(A)
        B = np.array(B)
    elif len(A) != len(B):
        raise Exception("your input vectors should have the same dimesion")
    return sum(list(A * B))

def mv(M,A):
    L = M * A
    l = len(A)
    e = np.ones(l)
    q = np.zeros(l)
    for i in range(l):
        q[i] = vv(L[i], e)
    return q

def mm(M, W):
    L,P = M.shape[1], W.shape[0]
    Q = np.zeros((L,P))
    for i in range(P):
        Q[i] = mv(M,W[i])
    return Q

def tm(T, M):
    L,P,P_ = T.shape[2], M.shape[0], M.shape[1]
    Q = np.zeros((L,P,P_))
    for i in range(L):
        Q[i] = mm(T[i], M)
    return Q



def E(G):
    e = np.zeros((4,4,4,4))
    e[0][1][2][3] = 1
    e[3][0][1][2] = 1
    e[2][3][0][1] = 1
    e[1][2][3][0] = 1
    e[3][2][1][0] = -1
    e[0][3][2][1] = -1
    e[1][0][3][2] = -1
    e[2][1][0][3] = -1
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    if i==j or i==k or  i==l or j==k or j==l or k==l:
                        e[i][j][k][l] = 0
    return np.sqrt(G)*e
print(E(1).shape)
def VV(A, B):
    if type(A) != np.ndarray:
        A = np.array(A)
    if len(A) != len(B):
        raise Exception('Incorrect shape')
    sum = 0
    for i in range(len(A)):
        sum += A[i] * B[i]
    return sum
#print(vecdot(A,B))

def MV(M, V):
    W =[]
    L = M.shape[0]
    for j in range(L):
        W.append(VV(M[j], V))
    return np.array(W)
#print(matvec(M, V))

def MM(E ,T):
    Q = np.zeros((E.shape[0],T.shape[1]))
    LL = E.shape[0]
    for p in range(LL):
        for q in range(LL):
            Q[:p] = MV(E[p],T[q])
    return np.array(Q)



def C_VMV(M, V1, V2):
    return VV(MV(M, V1), V2)


def C_TM(T, M):
    LLL = T.shape[0]
    I = np.zeros(T.shape[0])
    for i in range(T.shape[0]):
        I[i] = 1
    W = np.zeros(LLL)
    for i in range(LLL):
        for k in range(LLL):
            for j in range(LLL):
                W[i] = 0.5 * MM(T[k,j], M[j])
    return(np.array(W))


F = np.array([[0,3,5,7],
              [3,0,13,17],
              [5,-13,0,27],
              [-11,17,27,0]])


















def poly(z, y):
    x = np.array([1,2,3,4,5])
    def f(x, dtype=np.array):
        V = np.zeros((5,5))
        for i in range(5):
            for j in range(5):
                V[i][j] = x[i]**j
        return V

    a = sl.solve(f(x), y)

    def Z(z):
        q = list(np.zeros(5))
        for i in range(5):
            q[i] = z**(i)
        return np.array(q)

    return np.matmul(a, Z(z))


u = np.random.rand(5,5)
X = linspace(-500, 500, 1000)
figure()
for k in range(5):
    plot(X, poly(X, u[k]))
title("noice plots")

#--------------SAME FUNCTION, FEWER LINES-------------------

def PLY(z, y):
    x = np.array([1,2,3,4])
    V = np.vander(x)
    a = sl.solve(V, y)
    return np.polyval(a, z)

x = np.array([1,2,3,4,5])
V = np.vander(x)
A = V[:,1:5]
AA = sl.inv(np.matmul(A, A.T))
B = np.matmul(A.T, AA)
y = np.array([0.5, -2.0, 1.0, -0.5, 1.0])
c = np.matmul(B, y)

X = linspace(-500, 500, 1000)
figure()
plot(X, PLY(X, c))
title("noice plots")








u = np.random.rand(5,5)
X = linspace(-500, 500, 1000)
figure()
for k in range(5):
    plot(X, PLY(X, u[k]))
title("noice plots")

#--------------------


print(c)
print(B.shape)
print(B)
print(AA.shape)

print(AA)

print(V)
print(V[:,1:5])














def xi(u):
    i=0
    S = list(np.zeros(np.ndim(u)))
    while i < np.ndim(u):
        S.append((u[0] + u[i+1]+ u[i+2])/3)
        i += 1
    return np.array(S)
R = np.array([1,2,3])
print(xi(R))






















def co_corrector(A):
    if A.shape[0] != A.shape[1]:
        raise Exception("matrix has to be square")
    if type(A) != np.ndarray():
        A = np.array(A)
    V = A.T

    max_iter = 2
    step = 1.e-2
    tol = 0.3
    k = 1
    kk = 0
    L = A.shape[0]
    PA = []
    PV = []
    dphi_1 = np.zeros()
    dphi_2

    for i in range(L):
        for r in range(i+1,L):
            PA.append(np.arccos(np.matmul(A[i], A[r])/(np.abs(A[i])*np.abs(A[r]))))
            PV.append(np.arccos(np.matmul(V[i], V[r])/(np.abs(V[i])*np.abs(V[r]))))

    PA = np.array(PA)
    PV = np.array(PV)

    dphi_1 = np.zeros(PA.ndim)
    Theta_1 = np.zeros(PA.ndim)

    dphi_2 = np.zeros(PA.ndim)
    Theta_2 = np.zeros(PA.ndim)

    for i in range(PA.ndim-1):
        if abs(PA[i]) < tol and PA[i] > 0:
            dphi_1[i] = PA[i] + step
            while abs(dphi_1[i]) < tol:
                dphi_1[i] += step
            Theta_1[i] = dphi_1[i] - PA[i]


        elif abs(PA[i]) < tol and PA[i] < 0:
            dphi_1[i] = PA[i] - step
            while abs(dphi_1[i]) < tol:
               dphi_1[i] -=  step
            Theta_1[i] = dphi_1[i] - PA[i]


        elif abs(PA[i] - np.pi) < tol:
            dphi_1[i] = PA[i] - step
            while abs(PA[i] - np.pi) < tol:
                dphi_1[i] -= step
            Theta_1[i] = dphi_1[i] - PA[i]
pass



















#                if A[j] == A[i]:
#                    P.append(np.matmul(A[i + 1], A[j])/(np.abs(A[i+1])*np.abs(A[j])))

#                    else:
#                        np.matmul(A[i], A[j])/(np.abs(A[i])*np.abs(A[j]))















x = np.array([3,6,4,9,5,4])
y = np.array([409,6,9,4,5,7])
V = np.zeros((6,6))
for i in range(6):
    for j in range(6):
        V[i][j] = x[i] * y[j]
print(sl.det(V))


#-----------------------------------------------------

#Linear algebra methoda with SciPy and a bit about universal functions

#-----------------------------------------------------



A = np.array([[1,2,3],
              [4,5,6],
              [7,8,9]])
K = np.array([[1,2,8],
              [3,55,6],
              [7,8,89]])

(LU, piv) = sl.lu_factor(A)
x = np.zeros((3,3))
xi = sl.lu_solve((LU, piv), K)
print(xi)

print((LU, piv))
print("----------------")
print(LU)
print("----------------")
print(piv)

def heaviside(x):
    if x >= 0:
        return 1.
    else:
        return 0.
vheaviside = np.vectorize(heaviside)
print(vheaviside(np.array([-1, 2])))
print()

A = np.arange(13)
B= A[1:]
print(B)
a = B[:3]
b = B[3:6]
c = B[6:9]
d = B[9:12]
R = np.vstack([a,b,c,d])
print(R)
print(R[0,:])
print(R[2:])

x = np.array([5-k for k in range(5)])
print(x)
y = np.zeros((5,5))
for i in range(5):
    y[i] = x
print(y)

#class poly:
#    def __init__(self, V = None, b = None):
#        self.V = V
#        self.b = b
#        if V is None:
#            self.V = lambda x: self.V(x)
#        if b is None:
#            self.b = lambda x: self.b(x)

A = np.array([[1,2,3],
              [4,5,6],
              [7,8,9]])
a = np.array([1,2,3])

b = sl.solve(A,a)


#try:
#    def poly(y,z):
#        y = np.random.rand(5)
#        A =[]
#    x = np.array([5-k for k in range(5)])
#    V = np.zeros((5,5))
#    for i in range(5):
#       V[i] = x
#       V = np.array([[1,2,3,4,5     ],
#                      [6,7,8,9,10    ],
#                      [11,12,13,15,16],
#                      [17,18,19,20,21],
 #                     [22,23,24,25,26]])
#
#    a = sl.solve(V,y)
#    def Z(z):
##        q = np.zeros(5)
#        for j in range(5):
#            q[j] = z**j
#            return q
#        for k in range(5):
#            A.append(a[k]*Z(z)[k])
#            return sum(A)

#except:
#    Exception("some exception occured")
#else:
#    print("system nominal")

#------------------------------------------------------------------------------------------
#<-------------------------------------------------------------------------------------------
#def poly(z,y):
#    max_iter= 100
#    step = 0.5
#    STEP = np.zeros((6,6))
#    V = np.zeros((6,6))
#    k = 0
#    x = np.array([6-k for k in range(6)])
#    for i in range(6):
#        V[i] = x

#    for i in range(6):
#        for j in range(6):
#            STEP[i][j] = step

#    while k < max_iter:
 #       try:
 #           a = np.linalg.solve(V,y)
 #       except np.linalg.LinAlgError:
 #           STEP = np.zeros((6,6))
 #           W = V + STEP
 #           while np.linalg.det(W)==0:
 #               W =+ STEP
 #           a = np.linalg.solve(V, y)
 #   k =+ 1
 #   if k > max_iter and np.linalg.det(W) == 0:
 #       return "LOL, W is not yet nonsingular"
 #   else:
 #       q = np.zeros(6)
 #       for i in range(6):
 #           q[i] = z**i
 #       return np.matmul(a,q)
#u = np.array([1,2,3,7,4,6])
#print(poly(4,u))

def f(z, t):
    x = np.array([3,6,4,9,5,8])
    V = np.zeros((6,6))
    for i in range(6):
        V[i] = x
    l = sl.solve(V, t)
    q = np.zeros(6)
    for i in range(6):
        q[i] = z**i
    return np.matmul(q,l)

u = np.array([-2.0, 0.5, -2.0, 1.0, -0.5, 1.0])
print(f(2,u))


















#-------------------------------------------------------------------------------
#--------------------------------------------------------------------------------















#X = np.linspace(-5, 5, 1000)

#p = np.array([0.0, 0.5, 1.0, 1.5,2.0])
#figure()
#plot(X, poly(p,X), label= "p(z)")
#legend(loc = "lower center", fontsize ="small")













#x = np.linspace(-2*pi, 2*pi, 100)

#figure()
#plot(x, np.sin(x), label = "näsim sucks dick")
#plot(x, np.cos(x), label ="näsim  suger ännu större dick")
#legend(loc = "lower center", fontsize = "small")

#fig, ax = subplots()
#ax.plot(x, np.sin(x), label ="nsäim är en faggot")
#ax.plot(x, np.cos(x), label = "god does not exist")
#ax.legend(loc ="lower center", fontsize ="small")






















#def solvesys(z):
#    V = np.array([[1,2,3],
#              [4,5,6],
#              [7,8,9]])
#    b = np.array([7,8,9])
#
#    a = sl.solve(V, b)
#    if np.allclose(np.matmul(V, a), b)==False:
#        raise Exception("solution to system is wrong")

#    Z = np.array([z, z**2, z**3])
#    P = np.matmul(Z, a)
#    return P
#x = np.linspace(-2, 2, 1000)
#figure()
#plot(x, solvesys(x), "r")
#title("what")















































































































































