# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 08:52:31 2019

@author: User
"""

#Review of errors

#try:
#    #some code that might raise some erros>
#except ValueError:
#    print("oops, a ValueError occured)
#except TypeError:
#    print("we got a TypeError")
#except Exception:
#    print("some other kind of exception")

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np


#A = [1,2]
#B = [2,5]
#C = []
#l = len(A)
#for i in range(0,l):
#   C.append(A[i]*B[i])
#print(C)
#print(sum(C))



#def vecdot(A,B):
#    A = np.array([1., 2.])
#    B = np.array([2., 3.])
#    l = len(A)
#    C = []
#    for i in range(l):
#        C.append(A[i]*B[i])
#    print(C)
#    return sum(C)
#print(vecdot(A,B))


#def matvec(V, M):
#    V = np.zeros(2)
#    M = np.array([A,B])
#    MV = []
#    l = M.shape[0]
#    for i in range(l-1):
#        MV.append(vecdot(M[i],V))
#    return(matvec(V, M))
#print(MV)

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


class vectorspace:
	def __init__(self, basis):
		self.basis = np.array(basis)
#------------------------------------------------------
	def __add__(self, other):
		if type(other) == int:
            """
            det visar sig att denna funktion möjligtgör för en
            att lägga till en integer till varje  (elementwise) i basen förvektorrummet.
            SÅÅ detta betyder att dessa methods som hector gjorde
            skapade verkar på själva bas(vetorerna, matriserna,...etc)
            vi ser att additionen av en integer sker INUTI (self.basis + other)
            """
			return type(self)(self.basis + other)

		if type(other) == vectorspace:
			return type(self) (self.basis + other.basis)
#------------------------------------------------------
	def __radd__(self, other):
		if type(other) == int:
			return self + other



V = vectorspace([[1,2,4],[0,4,1],[2,3,9]])
W = vectorspace([[1,2,3], [5,6,8], [3,4,5]])
print((W+V + 5).basis)
print(V.__radd__(3))
print(vectorspace.__radd__(V,3).basis)



class Tensor:
    def __init__(self, basis, entry):
        if self.basis != np.ndarray:
            self.basis = np.array(self.basis)
        if self.entry != np.ndarray:
            self.entry = np.array(self.entry)
        if self.basis.shape != self.entry.shape:
            raise Exception('Incorrect shape of basis and tensor')


    def __add__(self, entry, other):
        if type(other) != np.ndarray:
            other = np.array(other)
        if other.shape != self.shape:
            raise Exception("do you want to add affinelly?")
        if type(other) == Tensor:
            return self.basis, type(self)(self.entry + other.entry)

    def


















































#def f(x,y,z):
#   if x < y < z:
#        return x**2 + y**2 - z**2
#    else:
#        print("condition not met")
#
#def g(x,y,z):
#    if f(x,y,z) == x+y+z - 100:
#        return (x, y, z)
#    else:
#        print("condition not met")
#c = 500
#for z in range(c):
#   for y in range(c):
#       for x in range(c):
#           if f(x,y,z) == 0 == x+y+z - 100:
#               print(x, ":", y, ":", z, ":", x*y*z)
#           else:
#               print("shit")
























#     if type(T) != np.ndarray:
#        T = np.array(T)
#     for s in range(LL):
#         if T[s].sahpe != M.shape:
#             raise Exception('Tens, Matr diff shape')


 #    Q = np.zeros(T.shape[0], T)


  #   for i in range(LL):
   #      for j in range(LL):
    #         for k in range(LL):
     #            Q[:i] = T[i][j][k]*M[j][k]
    # return np.array(Q)



#------------------------------------------------------------

W = np.array([[0., 3., 3.],
              [-3., 0., 3.],
              [-3., -3., 0.]])

e = np.array([[[0., 0., 0.],
              [-1., 0., 0. ],
               [0., -1., 0. ]],

             [[1., 0., 0.  ],
              [0., 0., 0. ],
              [0.,  0., -1.]],

             [[0., 1., 0.  ],
              [0., 0., 1],
              [0.,  0., 0.]]])

E = np.array([[0., 1.,  2.,  33.],
              [1., 5.,  6.,  66.],
              [2., 9.,  10., 84.],
              [3., 12., 13., 99.],
              [4., 15., 16., 46.],
              [5., 18., 19., 87.],
              [6., 21., 22., 45.],
              [7., 24., 25., 77.],
              [8., 27., 28., 44.]])

print(E.shape)
print(E[6:8, 1:3])
print("-")
print(E[1:8, 1:3])
print("-")

print(E[1:-1, 1:3] )
print("-")
print("-")
print("-")

print(E[3:7, 3:])
print("-")
print( E[6,: ])
print("plopp")
print(E[5:7, :])
print("new print")
print(E[1:3])

print("--------testing out the reshape function----------------")

v = np.array([k for k in range(16)])
M = v.reshape(4,-1)
print(M)

print("-------------testing out the transpose function vai reshape-------------")

N = M.T
print(N)

Q1 = [[1,2,3], [4,5,6]]
Q2 = np.array(Q1)
P = Q2.T
print(P)















































#def LTMV(L, v):
#    if type(L) != np.ndarray:
#        L = np.array(A)
#    if L.shape[0] != V.shape[0]:
#        raise Exception('Incorrect shape')
#    H = []
#   for k in range(0):
#        H.append(L[k][0] * V[k])
#        for s in range(1):
#            H.append(L[s][1] * V[s])
#            for t in range(2):
#                H.append(L[t][2] * V[t])
#    print(H)
#    return np.array(H)

#def LTMV(L, V):
#    if type(L) != np.ndarray:
#        L = np.array(V)
#    if L.shape[0] != V.shape[0]:
#        raise Exception('Incorrect shape')
#   H = []
#   for k in range(L.shape[0]):
#        for s in range(2):
#            H[:k] = L[s][k] * V[s]
#    print(H)
#    return np.array(H)


#Å = np.array([[1.,0.,0.],
 #            [1.,1.,0.],
  #           [1.,1.,1.]])

#Ä = np.array([1.,1.,1.])









#S = np.array([[3.,0.,0.],
 #            [0.,1.,0.],
  #           [0.,0.,1.]])

#v = np.array([8.,5.,6.])
#print(matvec(S,v))


#e = np.array([[ 0., 1., 1. ],
 #             [-1., 0., 1. ],
 #             [-1.,-1., 0.]])

#e = np.array([[[0., 1., 0.],
 #             [-1., 0., 0. ],
  #             [0.,  0., 0. ]],
#
 #            [[0., 0., 1.  ],
  #            [0., 0., 0. ],
   #           [-1.,  0., 0.]],
#
 #            [[0., 0., 0.  ],
  #            [0., 0., 1. ],
   #           [0.,  -1., 0.]]])
#def W(x):
#    def WW(x):
#        return np.array([[0, A_3(x), A_2(x)], [-A_3(x),0,A_1(x)], [-A_2(x),-A_1(x),0.]])
#    return WW(x) @ e
#print(W(3))











#Q = np.array([1.,2.,3.])
#P = np.array([3., 2., 1.])
#print(VV(Q, P))




































class CN:
    def __init__(self, Re, Im):
        if not isinstance(Re, int):
            raise TypeError("the real part should be a number")
        if not isinstance(Im, int):
            raise TypeError("the imaginary part should be a number")
        self.Re = Re
        self.Im = Im


    def __RET_Re__(self):
        return np.array([self.Re, self.Re])[0]
    def __RET_Im__(self):
        return np.array([self.Im, self.Im])[-1]


    def __CNadd__(self, other):
        if type(other) == type(self.Re):
            return np.array([(self.Re+other) , self.Im])

        if type(other) == type(self.Im):
            return np.array([self.Re , (self.Im+other)])

        if type(other) == CN:
            return np.array([self.Re + other.Re , self.Im + other.Im])


    def __CNmul__(self, other):
        a1, a2 = self.Re, self.Im
        b1, b2 = other.Re, other.Im
        return np.array([a1*b1 - a2*b2, a1*b2 + a2*b1])

    def __CNconj__(self):
        a1, a2 = self.Re, self.Im
        return np.array([a1, -a2])

    def __CNdiv__(self, other):
#        a1, a2 = self.Re, self.Im
#        b1, b2 = other.Re, other.Im
        return CN.__CNmul__(self,CN.__CNconj__(other))/abs(self)**2


#--------------------------------------CREATING PLOTS---------------------------------------------------

x1 = linspace(-2*pi, 2*pi, 100)
x2 = linspace(0, 2, 10)

# a figure is created
plot(x1, sin(x1), 'r')
plot(x1, cos(x1), 'b--')
title('My first figure')

# a new figure is created
figure()
plot(x2, sqrt(x2), 'go')
plot(x2, x2**2, 'y^')
title('My second figure')

#---------------------------------------Creating a splitscreen of two plots--------------------------

#the subplot command subplot(xyz)
#x = number of rows
#y =  number of kolums
#z = plot number
x1 = linspace(-2*pi, 2*pi, 100)
x2 = linspace(0, 2, 10)

 # not needed when using Spyder
figure()
subplot(211) # the first subplot using a 2 rows and 1 column grid
plot(x1, sin(x1), 'r')
plot(x1, cos(x1), 'b--')
title('My first subplot')

 # the second subplot using a 2 rows and 1 column grid
subplot(222)
plot(x2, sqrt(x2), 'go')
plot(x2, x2**2, 'y^')
title('My second subplot')















































def fac(x):
    if x == 0:
        return 1
    else:
        return x*fac(x-1)

def Re(x,n):
    return x**(2*n +1)/fac(2*n +1)
print(Re(1, 10))

for k in range(5):
    if Re(1,k) > 0.005:
        print("n", "=", k , ":", Re(1,k))
#------------------------------------------------TAYLOR POLYNOM FÖR ELEMTÄRA FUNKTIONER--------------------------------------------------

def VV(A, B):
    if type(A) != np.ndarray:
        A = np.array(A)
    if len(A) != len(B):
        raise Exception('Incorrect shape')
    sum = 0
    for i in range(len(A)):
        sum += A[i] * B[i]
    return sum
#------------------------------------------------------------------------e
def Taylor_e(x, n):
    def fac(s):
        if s == 0:
            return 1
        else:
            return s*fac(s-1)

    def VV(A, B):
        if type(A) != np.ndarray:
            A = np.array(A)
        if len(A) != len(B):
            raise Exception('Incorrect shape')
        sum = 0
        for i in range(len(A)):
            sum += A[i] * B[i]
        return sum

    M = np.zeros(n)
    N = np.zeros(n)
    for k in range(n):
        M[k] = x**(k)/fac(k)
    for j in range(n):
        N[j] = 1
    return VV(M,N)
print(Taylor_e(1,1000))
#------------------------------------------------------------------------cos
def Taylor_cos(x, n):
    def fac(s):
        if s == 0:
            return 1
        else:
            return s*fac(s-1)

    M = np.zeros(n)
    N = np.zeros(n)
    for k in range(n):
        M[k] = x**(2*k)/fac(2*k)*(-1)**(k)
    for j in range(n):
        N[j] = 1

    def VV(A, B):
        if type(A) != np.ndarray:
            A = np.array(A)
        if len(A) != len(B):
            raise Exception('Incorrect shape')
        sum = 0
        for i in range(len(A)):
            sum += A[i] * B[i]
        return sum
    return VV(M,N)

print(Taylor_cos(pi/3, 10))
#------------------------------------------------------------------------sin
def Taylor_sin(x,n):
    def fac(s):
        if s == 0:
            return 1
        else:
            return s*fac(s-1)

    M = np.zeros(n)
    N = np.zeros(n)
    for k in range(n):
        M[k] = x**(2*k +1)/fac(2*k +1)*(-1)**(k)
    for j in range(n):
        N[j] = 1

    def VV(A, B):
        if type(A) != np.ndarray:
            A = np.array(A)
        if len(A) != len(B):
            raise Exception('Incorrect shape')
        sum = 0
        for i in range(len(A)):
            sum += A[i] * B[i]
        return sum
    return VV(M,N)
print(Taylor_sin(pi/6,20))












