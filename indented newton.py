# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 13:18:46 2019

@author: User
"""
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
help(plt.contour)
#def fa(X):
#     x,y = X[0],X[1]
#     return x**2+y**2
#Q = np.array(optimize.fmin(fa,np.array([1,1]),retall=True)[-1])
#X = np.array([Q[-i,0] for i in range(1,10)])
#Y = np.array([Q[-i,1] for i in range(1,10)])
#fig, ax = plt.subplots()
#ax.plot(X,Y)
#print(Q)



def funcktion():



def contours(f, a,b,N,lineCount,x0):
     """
     Task 1 multivariable analysis
     """
     #a,b define your window in which contours are drawn

     Q = np.array(optimize.fmin(f,x0, retall=True)[-1])
     A = np.array([Q[-i,0] for i in range(Q.shape[0])])
     B = np.array([Q[-i,1] for i in range(Q.shape[0])])


     x = np.linspace(a,b, N)
     y = np.linspace(a,b, N)
     X,Y = np.meshgrid(x,y)
     Z = np.array([[f(np.array([X[i,j], Y[i,j]])) for i in range(X.shape[0])] for j in range(Y.shape[0])])

     fig, ax = plt.subplots()
     #plt.contour(X, Y, Z)
     ax.contour(X, Y, Z,lineCount)
     ax.plot(A,B)

     plt.title("contourplot")
     plt.xlabel("x")
     plt.ylabel("y")
     plt.show()


def fa(X):
     x,y = X[0],X[1]
     return x**2+y**2

print(contours(fa, -10,10,100, 20 ,np.array([7,3])))


import numpy as np
import matplotlib.pyplot as plot
import pylab

# List of points in x axis
XPoints     = []

# List of points in y axis
YPoints     = []

# X and Y points are from -6 to +6 varying in steps of 2
for val in range(-6, 8,2):
    XPoints.append(val)
    YPoints.append(val)

# Z values as a matrix
ZPoints = np.ndarray((7,7))

# Populate Z Values (a 7x7 matrix) - For a circle x^2+y^2=z
for x in range(0, len(XPoints)):
    for y in range(0, len(YPoints)):
        ZPoints[x][y] = (XPoints[x]* XPoints[x]**2) + (YPoints[y]*YPoints[y]**2)

# Print x,y and z values
print(XPoints)
print(YPoints)
print(ZPoints)

# Set the x axis and y axis limits
pylab.xlim([-10,10])
pylab.ylim([-10,10])

# Provide a title for the contour plot
plot.title('Contour plot')

# Set x axis label for the contour plot
plot.xlabel('X')

# Set y axis label for the contour plot
plot.ylabel('Y')

# Create contour lines or level curves using matplotlib.pyplot module
contours = plot.contour(XPoints, YPoints, ZPoints)

# Display z values on contour lines
plot.clabel(contours, inline=1, fontsize=10)

# Display the contour plot
plot.show()



#numbers = list(map(int, input('Enter a series of numbers: ').split(' ')))
#ind = [i for i, n in enumerate(numbers) if not n]
#print('Amount of numbers:', ind[0])
#

A = np.array([])
print(np.append(A, 8))




def number_of_numbers(list_):
     pre_zeros1 =np.array([])
     pre_zeros2 =np.array([])

     if isinstance(list_, list):
          def check(n):
               pre_zeros2 = np.append(pre_zeros1, n)
               return 1/n
          try:
               [check(item) for item in list_]
          except ZeroDivisionError:
               return len(pre_zeros2)-1, pre_zeros2
     else:
          raise Exception("your input must be a list containing numbers")

print(number_of_numbers([3,5,7,5,7,8,5,7,8,89,9,9,77,7,8,9,98,8,89,9,9,9,9,4,6,0,3,7,89,3,6,3,6,9,4,80]))


def number_of_numbers(list_):
     pre_zero =[]
     if isinstance(list_, list):
          def check(n):
               pre_zero.append(n)
               return 1/n
          try:
               [check(item) for item in list_]
          except ZeroDivisionError:
               return len(pre_zero)-1, pre_zero
     else:
          raise Exception("your input must be a list containing numbers")

print(number_of_numbers([3,5,7,5,7,8,5,7,8,89,9,9,77,7,8,9,98,8,89,9,9,9,9,4,6,0,3,7,89,3,6,3,6,9,4,80]))



def surface2D_curvature_tensor(J,ndiff, X):
     b = -np.matmul(J(*X).T,ndiff(*X))
     return b


def Test_F(x,y):

    return np.array([[1,0],[0,1],[3*x**3-3*x*y**2 -1, 3*x**2*y -y**3]])

def Test_n(x,y):
     return np.array([-(3*x**3-3*x*y**2 -1), -(3*x**2*y -y**3), 1])

def Test_ndiff(x,y):
     return np.array([[-3*x**2 + 3*y**2, 6*x*y], [-6*x**2, 3*y**2], [0,0.]])



b = surface2D_curvature_tensor(Test_F,Test_ndiff, np.array([-.01,-0.1]))
print(b, np.linalg.eigvals(b))











def ODE_solver(f, fp, t0, u0, te, N=None, h = 4):
    U = [[u0], [t0]]
    Y = [[u0], [t0]]
    delta_t = np.abs(te-t0)

    while np.abs(U[1][0] - U[1][-1]) <= delta_t:
        RadiusOfCurvature = fp(U[0][-1], U[1][-1])
        u_new = U[0][-1] + h*RadiusOfCurvature*f(U[0][-1], U[1][-1])
        t_new = U[1][-1] + h*RadiusOfCurvature

        U[0].append(u_new)
        U[1].append(t_new)

    while np.abs(Y[1][0] - Y[1][-1]) <= delta_t:
        u_new = Y[0][-1] + h*f(Y[0][-1], Y[1][-1])
        t_new = Y[1][-1] + h

        Y[0].append(u_new)
        Y[1].append(t_new)

    def y(x):
        return -3 -x - 2/(1+x) +5*np.exp(x)/(1+x)
    plt.figure()
    plt.plot(U[1], U[0], label = "curvature regulated euler",color =  "r")
    plt.plot(Y[1], Y[0], label = "standard euler",  color = "b")
    x = np.linspace(t0, te, (len(Y[1]) + len(U[1]))/2)
    plt.plot(x, y(x), label = "actual solution", color = "g")
    plt.title("comparing solutions to ODE")
    plt.legend(loc = "best")


def f(y,x):
    return y*(x/(1+x)) +1 +x

def fp(y,x):
    return  np.abs(x/(x + 1)*(y*(x/(1+x)) +1 +x) -x*y/(x + 1)**2 + y/(x + 1) + 1)/(1+ (y*(x/(1+x)) +1 +x)**2)

ODE_solver(f,fp, 0, 0, 10)



#
#
#
#def lvl(F, x0, y0, z0 h=1.e-20, error_message=True):
#     def Dx(F,a,b,h=h):
#          return (F(a+h,b) -F(a,b))/h
#     def Dy(F,a,b,h=h):
#          return (F(a,b+h) -F(a,b))/h
#     def yp(a,b):
#          return -Dx(F,a,b)/Dy(F,a,b)
#     def y_approx(t,a,b):
#          return b + yp(a,b)*(t-a)
#
#
#     xvals=[(x0,y0)]
#     k=0
#     while k < maxiter:
#          try:
#              new_yp = yp(*xvals[-1])
#
#          except ZeroDivisionError:
#               x_new +=
#
#
#
#     return None
#
#







def det_diag(matrix, show_details=False):
     a,b = matrix.shape[0], matrix.shape[1]

     def partial_column_zero(M, i,j):
          #(i: # of steps down, j: # steps right) within the matrix
          minor0 = M[i-1:a, j-1:b]
          minor1 = np.row_stack([minor0[0], np.array([minor0[0,0]*minor0[i] for i in range(1, minor0.shape[0])])])
          M[i-1:a, j-1:b] = np.row_stack([minor0[0], np.array([minor1[i] - minor0[0]*minor0[i,0] for i in range(1, minor1.shape[0])])])
          return M

     if show_details:
          for i in range(1,matrix.shape[0]):
               matrix = partial_column_zero(matrix, i,i)
          L = np.array([matrix[i,i] for i in range(matrix.shape[0])])
          return np.column_stack([matrix,L]), np.prod(L)
     else:
          if a != b:
               raise Exception("matrix must be square")
          Q = np.array([partial_column_zero(matrix, i,i)[i,i] for i in range(1,matrix.shape[0])])
          return np.prod(Q)

K = np.array([[0,1,0,0],
              [2,0,1,0],
              [0,2,0,1],
              [0,0,2,0]])
print(det_diag(K, show_details=True))















E = np.array([[3,8,4,6],
              [5,3,2,4],
              [7,1,4,3],
              [1,1,6,1]])









E = np.array([[3,8,4,6],
              [5,3,2,4],
              [7,1,4,3],
              [1,1,6,1]])
print(det_diag(E, show_details=True))
print(np.linalg.det(E))




#K = np.array([[1,4,5,8],
#              [1,1,9,2],
#              [1,7,1,6],
#              [6,2,1,1]])



def partial_column_zero(M, i,j):
     #(i: # of steps down, j: # steps right) within the matrix
     a,b = M.shape[0], M.shape[1]
     minor0 = M[i-1:a, j-1:b]
     minor1 = np.row_stack([minor0[0], np.array([minor0[0,0]*minor0[i] for i in range(1, minor0.shape[0])])])
     M[i-1:a, j-1:b] = np.row_stack([minor0[0], np.array([minor1[i] - minor0[0]*minor0[i,0] for i in range(1, minor1.shape[0])])])

     return M
K = np.array([[99,4,4.5,8],
              [1,1,-5,2],
              [1,7,1,88],
              [6,3,1,1]])

A = np.array([[3.,2,1,39,0],
             [2,3,1,34,0],
             [1,2,3,26,1]])
B = np.array([[3.,2,1],
             [2,3,1],
             [1,2,3]])
C = np.array([[1,2,-1],
             [2,-1,1],
             [-1,1,-2.]])
Q = np.array([[1,1,2],
              [1,-1,8],
              [2,3,1]])
W = np.array([[1,2,3.],
              [1,3,5],
              [1,4,6]])
print("......")
print(det_diag(C))
print(np.linalg.det(C))




E = np.column_stack([K, np.ones_like(K.T[0])])
print(K)
print("......")
print(K[0:4, 2:4])
print(K[2:4,2:4])

def subtraction(M):
     a,b  = M.shape[0], M.shape[1]

     M1 = M[1:a, 0:b]*M[0,0]
     M2 = np.array([M[0]*M[i,0] for i in range(1, a)])
     M[1:a, 0:b]  = M1 - M2
     return M





def PLOT(z0,N=None,pl=False, tol=1.e-10,max_iter=100000,eps1=0.1, eps2=1.e-16, matrix=None, normalised=False, quadraticform=False, epsilon_iter=False):
     """
     task 2
     """
     def convergence_tests1(z0, tol=tol, max_iter=max_iter, matrix=matrix, normalised=normalised, quadraticform=quadraticform):
          if not matrix:
               K = np.array([[1,1,2.],[-3,4,3],[2,3,1]])
          else:
               K = matrix
          def Qform(M, v):
               vecmat = v*v.reshape(-1,1)
               return np.trace(np.matmul(M, vecmat))

          zmax  = 1/tol
          k     = 0
          zvals = [z0]
          znorm = np.linalg.norm(zvals[-1])
          qvals = [Qform(K,z0)]

          while k < max_iter and znorm < zmax:
               new = np.matmul(K, zvals[-1])

               if normalised:
                    new = new/(np.linalg.norm(new))
                    zvals.append(new)

                    if quadraticform == False:
                         delta = np.linalg.norm(zvals[-1] - zvals[-2])
                         if delta < tol:
                              return "z convergent, delta = {delta}", zvals[-1], k
                    else:
                         qnew = Qform(K, new)
                         qvals.append(qnew)
                         deltaq = np.linalg.norm(qvals[-1] - qvals[-2])
                         if deltaq < tol:
                              return f"q convergent, deltaq = {deltaq}", qvals[-1], k
               else:
                    zvals.append(new)
                    delta = np.linalg.norm(zvals[-1] - zvals[-2])
                    if delta < tol:
                         return f"z convergent delta = {delta}", np.linalg.norm(zvals[-1]),k
               k +=1

          if quadraticform == False:
               return f"z divergent", zvals[-1], k
          else:
               return "q divergent", qvals[-1], k

     if pl:
          if not N and normalised == False == quadraticform:
               raise Exception("you must input N and turn on normalised- and quadraticform key words")
          x_axis = np.linspace(eps1,eps2,N)
          if epsilon_iter:
               yznorm_list = np.array([convergence_tests1(z0, normalised=True,tol = i)[-1] for i in x_axis])
               yzqnorm_list = np.array([convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-1] for i in x_axis])
               difference = np.array([np.abs(convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-1]- convergence_tests1(z0, normalised=True,tol = i)[-1]) for i in x_axis])
          else:
               yznorm_list = np.array([convergence_tests1(z0, normalised=True,tol = i)[-2] for i in x_axis])
               yzqnorm_list = np.array([convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-2] for i in x_axis])

          fig, ax = plt.subplots()
          ax.plot(x_axis, yznorm_list)
          plt.title("normalised z, iterations vs epsilon")
          plt.xlabel("epsilon, interval = [10^(-16),0.1]")
          plt.ylabel("iterations")

          fig, ax = plt.subplots()
          ax.plot(x_axis, yzqnorm_list)
          plt.title("normalised z with q, iterations vs epsilon")
          plt.xlabel("epsilon, interval = [10^(-16),0.1]")
          plt.ylabel("iterations")

          if epsilon_iter:
               fig, ax = plt.subplots()
               ax.plot(x_axis, difference)
               plt.title("difference in iterations")
               plt.xlabel("epsilon, interval = [10^(-16),0.1]")
               plt.ylabel("iterations")

          plt.show()
     else:
          return convergence_tests1(z0)



print(PLOT(np.array([8,3,13]), N =1000, normalised=False,  epsilon_iter=True))





def lin_dep(Set, margin=1.e-5, error_notification="raise_excpetion"):
    """
    This function checks if a set of vectors are linearly dependent or not.

    One of the conditions that makes a set of vectors linearly independent is that Xa = 0
    results in a vector a whose components all are zero.
    if they would have been anything else then they would be linearly dependent

    On INPUTS:
        Set (list):
            a list of arrays
        error_notification (string/bool):
            can either be the string "raise_excpetion" or the boolena  value None
    ON RETURN:


    """
    Dim = []
    for i, s in enumerate(Set):
        if type(s) != np.ndarray or np.ndim(s) != 1:
            if not error_notification:
                Set = Set[:i] + Set[-(len(Set) -(i+1)):]
            else:
                raise Exception("Set must only contain vectors")

    for s in Set:
        Dim.append(s.shape[0])

    if max(Dim) == min(Dim) < len(Dim):
        return "linearly dependent 0: amount of vectors is larger than their dimension"
    else:
        m_Dim = max(Dim)
        ind = Dim.index(m_Dim)
        M = np.zeros((len(Set), Set[ind].shape[0]))
        for i, s in enumerate(Set):
            s_z = np.zeros(m_Dim - Dim[i])
            s = np.hstack([s, s_z])
            M[i] = s

        zero = np.zeros(m_Dim)
        try:
            conclusion = np.linalg.solve(M, zero)
        except sl.LinAlgError:
            return np.linalg.det(M),"linearly dependent 1", True

        if np.linalg.norm(conclusion) < margin:
            return f"det = {np.linalg.det(M)}",conclusion, "linearly independent 2", False
        else:
            return M,conclusion, "linearly dependent 3", True


K = np.array([[1,0,0,0],
              [1,-1,0,2],
              [1,0,-1,0],
              [0,2,-1,1]])
Set = list(K.T)
print(lin_dep(Set))



def x(t):
     return np.sin(2*t)
def y(t):
     return 0.5*(np.sin(t))**2

T = np.linspace(0,10, 100)
X = [x(t) for t in T]
Y = [y(t) for t in T]
fig, ax = plt.subplots()
ax.plot(X, Y)
plt.show()



print(1/np.sqrt(12))
print(0.28867513*3)
c = 12**2+ 2*2+ 24**2+ 35**2

print(c*24)









#
#def least_square(transformation_matrix, y):
#     """
#     the kolumns of the transformation_matrix are the basis for its outputspace
#     """
#     def m(a,b):
#               return np.matmul(a,b)
#     def n(a):
#          return np.linalg.norm(a)
#
#     def scalar_gen(y,basisvector):
#          t = m(y, basisvector)/n(basisvector)**2
#          return t
#
#     Mlist = list(transformation_matrix.T)
#     scalars = np.array([scalar_gen(y,basisvector) for basisvector in Mlist])
#     k = m(transformation_matrix, scalars) #k is a scaled version of y_approx
#     xi = scalar_gen(y, k)
#     y_approx = xi*k
#
#     x_min = np.linalg.solve(A, y_approx)
#     return x_min
#
print(1/np.sqrt(3.333333333333333333333333333))

def LS(A, y):
     def m(a,b):
               return np.matmul(a,b)
     def n(a):
          return np.linalg.norm(a)
     def ON_basis_gen(matrix, orthogonality_check=False, automatic_check=False, error_tol=1.e-10, QR_factorisation=False):
          veclist,newbasis = list(matrix), []
          def orth_check(Matrix):

               list_= list(Matrix)
               dot_matrix = np.array([[m(item1, item2) for item1 in list_] for item2 in list_])
               if (dot_matrix - np.eye(dot_matrix.shape[0]) < error_tol).all():
                    return True
               else:
                    error = dot_matrix - np.eye(dot_matrix.shape[0])
                    return False, np.max(error), np.min(error)
          def motor(vector, ind):
               if ind == 0:
                    new = vector/n(vector)
               else:
                    L = np.array([newbasis[i]*m(newbasis[i],vector) for i in range(len(newbasis))])
                    projections = np.sum(L, axis=0)
                    NEW = vector - projections
                    new = NEW/n(NEW)
               newbasis.append(new)
          [motor(vector, ind) for ind,vector in enumerate(veclist)]
          newbasis_matrix = np.array(newbasis)
          if QR_factorisation:
               Q, Qt, A = newbasis_matrix.T,newbasis_matrix, matrix.T
               R = m(Qt, A)
               QRA = np.array([Q,R,A])
               return QRA
          else:
               if orthogonality_check and automatic_check == False:
                    return orth_check(newbasis_matrix)
               elif automatic_check:
                    return orth_check(matrix)
               else:
                    return newbasis_matrix
     QRA = ON_basis_gen(A, QR_factorisation=True)
     try:
          Qty = m(QRA[0].T, y)
     except ValueError:
          QRA = ON_basis_gen(A.T, QR_factorisation=True)
          Qty = m(QRA[0].T, y)
     x = np.linalg.solve(QRA[1], Qty)
     return x

A = np.array([[1,1,1],
              [1,2,2],
              [2,1,2],
              [2,1,0.]])
y = np.array([1,0,1,0.])
print(LS(A,y))
A = np.array([[1,1.],
              [4,-1],
              [3,2]])

y = np.array([6,8,5.])
print(LS(A, y))




A = np.array([[1,1,2.],
              [1,2,1],
              [2,1,1],
              [2,2,1]])

y = np.array([1,-1,1,-1])
print(LS(A.T,y,QR_method=True, normal_method=True))
x = np.linalg.lstsq(A,y,rcond=None)
print(x[0])




















#def qr_null(A, tol=None):
#    Q, R, P = qr(A.T, mode='full', pivoting=True)
#    tol = np.finfo(R.dtype).eps if tol is None else tol
#    rnk = min(A.shape) - np.abs(np.diag(R))[::-1].searchsorted(tol)
#    return Q[:, rnk:].conj()
#
#A = np.array([[3,2,5,-1.],
#              [1,-2,3,6],
#              [2,-4,6,12]])
#Z = qr_null(A)
#print(Z)
#print(np.matmul(A, Z))







def Grahm_Schmidt(matrix, lin_dep_inspection=False):
     def lin_dep(Set, margin=1.e-5, error_notification="raise_excpetion"):
         """
         This function checks if a set of vectors are linearly dependent or not.

         One of the conditions that makes a set of vectors linearly independent is that Xa = 0
         results in a vector "a" whose components all are zero.
         if only one of the components would have been anything else then the set of vectors
         would be linearly dependent

         On INPUTS:
             Set (list):
                 a list of arrays

             error_notification (string/bool):
                 can either be the string "raise_excpetion" or the boolena  value None
         ON RETURN:
              lol
         """
         if len(Set) == 2:
              if np.abs(np.matmul(Set[0], Set[1])) < margin:
                   return "linearly dependent 0", True
              else:
                   return "linearly independent 0", False
         Dim = []
         for i, s in enumerate(Set):
             if type(s) != np.ndarray or np.ndim(s) != 1:
                 if not error_notification:
                     Set = Set[:i] + Set[-(len(Set) -(i+1)):]
                 else:
                     raise Exception("Set must only contain vectors")

         for s in Set:
             Dim.append(s.shape[0])

         if max(Dim) == min(Dim) < len(Dim):
             return "linearly dependent 0: amount of vectors is larger than their dimension",True
         else:
             m_Dim = max(Dim)
             ind = Dim.index(m_Dim)
             M = np.zeros((len(Set), Set[ind].shape[0]))
             for i, s in enumerate(Set):
                 s_z = np.zeros(m_Dim - Dim[i])
                 s = np.hstack([s, s_z])
                 M[i] = s

             zero = np.zeros(m_Dim)
             try:
                 conclusion = np.linalg.solve(M, zero)
             except sl.LinAlgError:
                 return np.linalg.det(M),"linearly dependent 1", True

             if np.linalg.norm(conclusion) < margin:
                 return f"det = {np.linalg.det(M)}",conclusion, "linearly independent 2", False
             else:
                 return M,conclusion, "linearly dependent 3", True

     veclist = list(matrix)

     if lin_dep_inspection:
          lindep_check = lin_dep(veclist)
          if lindep_check[-1]:
               raise Exception("the set of vectors are linearly dependent")
     else:
          newbasis = []
          def m(a,b):
                    return np.matmul(a,b)
          def n(a):
               return np.linalg.norm(a)

          for ind, vector in enumerate(veclist):
               if ind == 0:
                    new = vector/n(vector)
               else:
                    projections = 0
                    for i in range(len(newbasis)):
                         projections += newbasis[i]*m(newbasis[i],vector)
                    NEW = vector - projections
                    new = NEW/n(NEW)
               newbasis.append(new)

     newbasis_matrix = np.array([i for i in newbasis])
     return newbasis_matrix




def ON_basis_gen(matrix, orthogonality_check=False, automatic_check=False, error_tol=1.e-10, QR_factorisation=False):
     """
     matrix is a matrix whose rows are linearly independent
     vectors to be turned into an ON-basis
     """
     veclist = list(matrix)
     newbasis = []
     def orth_check(Matrix):
          """
          This fucntion check for the pairwise orthogonality of the new basis
          """
          list_ = list(Matrix)
          dot_matrix = np.array([[m(item1, item2) for item1 in list_] for item2 in list_])
          if (dot_matrix - np.eye(dot_matrix.shape[0]) < error_tol).all():
               return True
          else:
               error = dot_matrix - np.eye(dot_matrix.shape[0])
               return False, np.max(error), np.min(error)

     def m(a,b):
               return np.matmul(a,b)
     def n(a):
          return np.linalg.norm(a)

     def motor(vector, ind):
          if ind == 0:
               new = vector/n(vector)
          else:
               L = np.array([newbasis[i]*m(newbasis[i],vector) for i in range(len(newbasis))])
               projections = np.sum(L, axis=0)
               NEW = vector - projections
               new = NEW/n(NEW)
          newbasis.append(new)

     [motor(vector, ind) for ind,vector in enumerate(veclist)]
     newbasis_matrix = np.array(newbasis)

     if QR_factorisation:
          Q, Qt, A = newbasis_matrix.T, newbasis, matrix.T
          R = m(Qt, A)
          QRA = np.array([Q,R, A])
          print("Q,R,A =")
          return QRA
     else:
          if orthogonality_check and automatic_check == False:
               return orth_check(newbasis_matrix)

          elif automatic_check:
               return orth_check(matrix)
          else:
               return newbasis_matrix

U = np.array([[1,1,1,1.],
              [1,1,1,0]])
print(ON_basis_gen(U))


K = np.array([[1,0,0,0],
              [1,-1,0,2],
              [1,0,-1,0],
              [0,2,-1,1]])
Kp = ON_basis_gen(K)
pK = np.linalg.inv(Kp)
print(Kp)
print(np.matmul(Kp, np.array([1,2,3,4.])))

print(ON_basis_gen(Kp, automatic_check=True))
print(Kp)


K = np.array([[1,1,1,1],
              [2.0001,2,2,2],
              [2,3,1,6]])
Kp = ON_basis_gen(K, orthogonality_check=True)
print(Kp, np.matmul(Kp, np.array([1,2,3,4.])))
print(np.matmul(Kp[0], Kp[2]))

A = np.array([[3,2,5,-1.],
              [1,-2,3,6],
              [2,3,-4,-3],
              [4,-1,2,1]])
C = Grahm_Schmidt(A)
print(C)

print(np.matmul(C[1], C[0]))






















fig = plt.figure(figsize = (4,2))
sld_ax = plt.axes([0.2, 0.9, 0.6, 0.1])# axes for slider
ax = plt.axes([0.1, 0.15, 0.8, 0.7])     # axes for the plot
sld = mpw.Slider(sld_ax, 'amp', 0., 5.)
x = np.linspace(-2*np.pi, 2*np.pi, 200)
ax.set_ylim(-5.5, 5.5)
lines, = ax.plot(x, sld.val*np.sin(x))
def update_amplitude(val):
     lines.set_ydata(val*np.sin(x))
sld.on_changed(update_amplitude)




fig = plt.figure()
slda_ax = plt.axes([0.05, 0.5, 0.3, 0.03])
# axes for a
sldd_ax =plt.axes([0.05, 0.4, 0.3, 0.03]) # axes for d
ax = plt.axes([0.5, 0.1, 0.4, 0.8])        # axes for plot
t = np.linspace(0, 14*np.pi, 400)
a0, d0 = 1., 0.05
slda = mpw.Slider(slda_ax, 'a', 0, 2, valinit=a0, valfmt="%1.2f")
sldd = mpw.Slider(sldd_ax, 'd', 0, 0.1, valinit=d0,  valfmt="%1.3f")
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_aspect(1.0)
# the aspect ratio should be 1 here
line = ax.plot(np.cos(a0*t)*np.exp(-d0*t), np.sin(t)*np.exp(-d0*t))
def update_curve(val):
     a = slda.val
     d = sldd.val
     line.set_xdata(np.cos(a*t)*np.exp(-d*t))
     line.set_ydata(np.sin(t)*np.exp(-d*t))
     slda.on_changed(update_curve)
     sldd.on_changed(update_curve)








def PLOT(z0,N=None,pl=False, tol=1.e-10,max_iter=100000,eps1=0.1, eps2=1.e-16, matrix=None, normalised=False, quadraticform=False, epsilon_iter=False):
     """
     task 2
     """
     def convergence_tests1(z0, tol=tol, max_iter=max_iter, matrix=matrix, normalised=normalised, quadraticform=quadraticform):
          if not matrix:
               K = np.array([[1,1,2.],[-3,4,3],[2,3,1]])
          else:
               K = matrix
          def Qform(M, v):
               vecmat = v*v.reshape(-1,1)
               return np.trace(np.matmul(M, vecmat))

          zmax  = 1/tol
          k     = 0
          zvals = [z0]
          znorm = np.linalg.norm(zvals[-1])
          qvals = [Qform(K,z0)]

          while k < max_iter and znorm < zmax:
               new = np.matmul(K, zvals[-1])

               if normalised:
                    new = new/(np.linalg.norm(new))
                    zvals.append(new)

                    if quadraticform == False:
                         delta = np.linalg.norm(zvals[-1] - zvals[-2])
                         if delta < tol:
                              return "z convergent, delta = {delta}", zvals[-1], k
                    else:
                         qnew = Qform(K, new)
                         qvals.append(qnew)
                         deltaq = np.linalg.norm(qvals[-1] - qvals[-2])
                         if deltaq < tol:
                              return f"q convergent, deltaq = {deltaq}", qvals[-1], k
               else:
                    zvals.append(new)
                    delta = np.linalg.norm(zvals[-1] - zvals[-2])
                    if delta < tol:
                         return f"z convergent delta = {delta}", np.linalg.norm(zvals[-1]),k
               k +=1

          if quadraticform == False:
               return f"z divergent", zvals[-1], k
          else:
               return "q divergent", qvals[-1], k

     if pl:
          if not N and normalised == False == quadraticform:
               raise Exception("you must input N and turn on normalised- and quadraticform key words")
          x_axis = np.linspace(eps1,eps2,N)
          if epsilon_iter:
               yznorm_list = np.array([convergence_tests1(z0, normalised=True,tol = i)[-1] for i in x_axis])
               yzqnorm_list = np.array([convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-1] for i in x_axis])
               difference = np.array([np.abs(convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-1]- convergence_tests1(z0, normalised=True,tol = i)[-1]) for i in x_axis])
          else:
               yznorm_list = np.array([convergence_tests1(z0, normalised=True,tol = i)[-2] for i in x_axis])
               yzqnorm_list = np.array([convergence_tests1(z0, normalised=True, quadraticform=True,tol = i)[-2] for i in x_axis])

          fig, ax = plt.subplots()
          ax.plot(x_axis, yznorm_list)
          plt.title("normalised z, iterations vs epsilon")
          plt.xlabel("epsilon, interval = [10^(-16),0.1]")
          plt.ylabel("iterations")

          fig, ax = plt.subplots()
          ax.plot(x_axis, yzqnorm_list)
          plt.title("normalised z with q, iterations vs epsilon")
          plt.xlabel("epsilon, interval = [10^(-16),0.1]")
          plt.ylabel("iterations")

          if epsilon_iter:
               fig, ax = plt.subplots()
               ax.plot(x_axis, difference)
               plt.title("difference in iterations")
               plt.xlabel("epsilon, interval = [10^(-16),0.1]")
               plt.ylabel("iterations")

          plt.show()
     else:
          return convergence_tests1(z0)



print(PLOT(np.array([8,3,13]), N =1000,pl=True, normalised=False,  epsilon_iter=True))

















P = [np.array([1,1.])]

def append_diff(Vector, list_):
     list_.append(Vector)
     delta = np.linalg.norm(list_[-1] - list_[-2])
     return delta
print(append_diff(np.array([2,2.]),P))
print(P)

def PLOT(z0,N=None,pl=False, tol=1.e-4,max_iter=1000,eps1=0.1, eps2=1.e-16, matrix=None, normalised=False, quadraticform=False, epsilon_iter=False):
     def convergence_tests(z0, tol=tol, max_iter=max_iter, matrix=matrix, normalised=normalised, quadraticform=quadraticform, epsilon_iter=epsilon_iter):
          if not matrix:
               K = np.array([[1,1,2.],[-3,4,3],[2,3,1]])
          else:
               K = matrix
          def Qform(M, v):
               vecmat = v*v.reshape(-1,1)
               return np.trace(np.matmul(M, vecmat))

          zmax  = 1/tol
          k     = 0
          zvals = [z0]
          znorm = np.linalg.norm(zvals[-1])
          qvals = [Qform(K,z0)]

          def Norm(v):
               abs_v = np.linalg.norm(v)
               w = v/abs_v
               return w
          def append_diff(Vector, list_, tol=tol, epsilon_iter=epsilon_iter):
               list_.append(Vector)

               if isinstance(Vector,np.ndarray):
                    delta = np.linalg.norm(list_[-1] - list_[-2])
               elif isinstance(Vector,float) or isinstance(Vector,int):
                    delta = np.abs(list_[-1] - list_[-2])

               if delta < tol:
                    if epsilon_iter:
                         return k
                    else:
                         return "convergent", list_[-1]

          while k < max_iter and znorm < zmax:
               new = np.matmul(K,zvals[-1])
               if normalised and quadraticform:
                    new = Norm(new)
                    qnew = Qform(K, new)
                    append_diff(new, zvals)
                    append_diff(qnew, qvals)
               elif normalised and quadraticform ==False:
                    new = Norm(new)
                    append_diff(new, zvals)
               elif normalised == False and quadraticform==True:
                    raise Exception("qform functionality requirs setting normalised = True" )
               else:
                    append_diff(new, zvals)
               k +=1
               print(zvals[-1])
          return f"divergent, last value computed = {zvals[-1]}"

     if pl:
          if not N:
               raise Exception("you must input N")
          x_axis = np.linspace(eps1,eps2,N)
          yz_list = np.array([convergence_tests(z0, epsilon_iter=True, tol = i) for i in x_axis])
          yznorm_list = np.array([convergence_tests(z0, epsilon_iter=True, normalised=True,tol = i) for i in x_axis])
          yzqnorm_list = np.array([convergence_tests(z0, epsilon_iter=True, normalised=True,tol = i) for i in x_axis])

          fig, ax = plt.subplots()
          ax.plot(x_axis, yz_list, label ="iterations vs epsilon")
          ax.plot(x_axis, yznorm_list, label ="normalised z, iterations vs epsilon")
          ax.plot(x_axis, yzqnorm_list, label ="normalised z with q, iterations vs epsilon")

          plt.show()
     else:
          return convergence_tests(z0)

print(PLOT(np.array([8,3,12]), normalised=True, quadraticform=True))









A = np.array([[1,1,2.],
              [1,2,1],
              [2,1,1],
              [2,2,1]])
y = np.array([1,-1,1,-1])
x = np.linalg.lstsq(A,y,rcond=None)
print(x[0])


def plott(val1,val2, N, plottt=None, a=None, method="numpy", method_difference=False, matrix=None):
     """
     Task 1
     This function can calculate the least square with respect to a matrix K, and a vector y.
     the plotting function is built inside the fcuntion and you can choose which of the built in numpy
     or scipy least square to use when calculating the minimum, but it's set by default as method = "numpy"
     """

     def LS(A, y):
          def m(a,b):
                    return np.matmul(a,b)
          def n(a):
               return np.linalg.norm(a)
          def ON_basis_gen(matrix, orthogonality_check=False, automatic_check=False, error_tol=1.e-10, QR_factorisation=False):
               veclist,newbasis = list(matrix), []
               def orth_check(Matrix):

                    list_= list(Matrix)
                    dot_matrix = np.array([[m(item1, item2) for item1 in list_] for item2 in list_])
                    if (dot_matrix - np.eye(dot_matrix.shape[0]) < error_tol).all():
                         return True
                    else:
                         error = dot_matrix - np.eye(dot_matrix.shape[0])
                         return False, np.max(error), np.min(error)
               def motor(vector, ind):
                    if ind == 0:
                         new = vector/n(vector)
                    else:
                         L = np.array([newbasis[i]*m(newbasis[i],vector) for i in range(len(newbasis))])
                         projections = np.sum(L, axis=0)
                         NEW = vector - projections
                         new = NEW/n(NEW)
                    newbasis.append(new)
               [motor(vector, ind) for ind,vector in enumerate(veclist)]
               newbasis_matrix = np.array(newbasis)
               if QR_factorisation:
                    Q, Qt, A = newbasis_matrix.T,newbasis_matrix, matrix.T
                    R = m(Qt, A)
                    QRA = np.array([Q,R,A])
                    return QRA
               else:
                    if orthogonality_check and automatic_check == False:
                         return orth_check(newbasis_matrix)
                    elif automatic_check:
                         return orth_check(matrix)
                    else:
                         return newbasis_matrix
          QRA = ON_basis_gen(A, QR_factorisation=True)
          try:
               Qty = m(QRA[0].T, y)
          except ValueError:
               QRA = ON_basis_gen(A.T, QR_factorisation=True)
               Qty = m(QRA[0].T, y)
          x = np.linalg.solve(QRA[1], Qty)
          return x

     def minimiser(a, method, matrix=matrix):
          def F(x):
               sigma = np.matmul(K,x)
               delta = sigma-y
               Norm = np.linalg.norm(delta)
               return Norm

          if not matrix:
               K = np.array([[1,1,2.],[1,2,1],[2,1,1],[2,2,1]])
          else:
               K = matrix

          if not a:
               y = np.array([1,-1,1,-1])
          else:
               y = np.array([1,a,1,a])

          if method == "numpy":
               X = np.linalg.lstsq(A,y,rcond=None)[0]
          elif method == "scipy":
               X = optimize.fmin(F, x0=np.zeros(3))
          elif method == "LS":
               X = LS(A, y)

          return np.linalg.norm(X)

     if not plottt and a != None:
          return minimiser(a)
     elif plottt != None:
          interval = np.linspace(val1,val2,N)
          fig, ax = plt.subplots()

          if method_difference:
               Y = np.array([np.abs(minimiser(x, "numpy") - minimiser(x, "scipy")) for x in list(interval)])
          else:
               Y = np.array([minimiser(x, method) for x in list(interval)])
     ax.plot(interval, Y, label ="norm vsersus a")
     plt.show()




plott(-1,5,100, plott=True, method="numpy")

plott(-1,5, 100, plott=True, method_difference=True)











fig = plt.figure(figsize = (4,2))
sld_ax = plt.axes([0.2, 0.9, 0.6, 0.1])# axes for slider
ax = plt.axes([0.1, 0.15, 0.8, 0.7])     # axes for the plot
sld = mpw.Slider(sld_ax, 'amp', 0., 5.)
x = np.linspace(-2*np.pi, 2*np.pi, 200)
ax.set_ylim(-5.5, 5.5)
lines, = ax.plot(x, sld.val*np.sin(x))
def update_amplitude(val):
     lines.set_ydata(val*np.sin(x))
sld.on_changed(update_amplitude)


fig = plt.figure()
slda_ax = plt.axes([0.05, 0.5, 0.3, 0.03])
# axes for a
sldd_ax =plt.axes([0.05, 0.4, 0.3, 0.03]) # axes for d
ax = plt.axes([0.5, 0.1, 0.4, 0.8])        # axes for plot
t = np.linspace(0, 14*np.pi, 400)
a0, d0 = 1., 0.05
slda = mpw.Slider(slda_ax, 'a', 0, 2, valinit=a0, valfmt="%1.2f")
sldd = mpw.Slider(sldd_ax, 'd', 0, 0.1, valinit=d0,  valfmt="%1.3f")
ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
ax.set_aspect(1.0)
# the aspect ratio should be 1 here
line, = ax.plot(np.cos(a0*t)*np.exp(-d0*t), np.sin(t)*np.exp(-d0*t))
def update_curve(val):
     a = slda.val
     d = sldd.val
     line.set_xdata(np.cos(a*t)*np.exp(-d*t))
     line.set_ydata(np.sin(t)*np.exp(-d*t))
     slda.on_changed(update_curve)
     sldd.on_changed(update_curve)



def p(r,v):
     def x(r,v):
          return r*np.cos(v)
     def y(r,v):
          return r*np.sin(v)
     X = x(r,v)
     Y = y(r,v)
     if X**2 + Y**2 > 1:
          return  np.log(X**2 + Y**2 -1)
     else:
          pass



#
#x = np.linspace(-1,1,100)
#y = np.linspace(-1,1,100)
#X, Y = np.meshgrid(x,y)
#L = np.array([[np.array([X[i,j], Y[i,j]]) for i in range(X.shape[0])] for j in range(Y.shape[0])])
#A = np.array([[p(*L[i,j]) for i in range(100)] for j in range(100)])
#
#fig = plt.figure()
#ax = fig.gca(projection ="3d")
#ax.plot_surface(X,Y,A,alpha=0.9)
#
##ax.plot_wireframe(X,Y,A, rstride=5, cstride=5)
#
##plot contour projection of each axis face
#ax.contour(X,Y, A, zdir="x", offset = -1)
#ax.contour(X,Y, A, zdir="y", offset = 1)
#ax.contour(X,Y, A, zdir="z", offset = -2)
#
##set limits
#ax.set_xlim3d(-1,1)
#ax.set_ylim3d(-1,1)
#ax.set_zlim3d(-2,2)
#
##set labels
#ax.set_xlabel("X axis")
#ax.set_ylabel("Y axis")
#ax.set_zlabel("Z axis")
#
#
#
#
#
#
#
#
#
#
##X,Y,Z = axes3d.get_test_data(0.05)
##X = np.linspace(-1,1,100)
##Y = np.linspace(-1,1,100)
#fig = plt.figure()
#ax = fig.gca(projection="3d")
#ax.plot_surface(X,Y,A,alpha=0.9)
#
#
#





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


def p1(x,y):
     return x**2 + y**2
def p1x(x,y):
     return 2*x
def p1y(x,y):
     return 2*y
def p2(x,y):
     return x+y +4
def p3(x,y):
     return 2*x + 5*y + 10

def level_sets(f, x1,x2, y1,y2, N, fx=None, fy=None, g=None, lvl_vals=None, function_list=None, tol=1.e-3, plot_style="3D", plane_face="xyz", appearence="surface", wireframedensity=2, coloroption="plasma", opacity=0.9):
     """
     Ok, listen up: this function creates level sets of a function f:R^2 --> R
     using meshgrid.
     PARAMETERS:
          f (function):
               this is a function, (mathematical in nature) of two variables
          lvl_vals (list):
               this is a list containing you surfaces of constant z-value
          x1,x2,y1,y2 (integers):
               this is the rectangle in which the function creates level sets
          N (integer):
               this is the resolution of the boundaries of the meshgrid
          tol (float):
               this value lets the function have some room when deciding if the function
               f equals any given level-value. two floats are almost never equal when expressed
               numerically in python so we must give the function some margin for error.
          plot_style (string):
               For 3D plotting:
                    plot_style = "3D"
                         appearance = "surface":
                              creates the literal surface of the function
                         plane_face:
                              creates level sets on the planes given by "plane_face":
                              plane_face can have the values:
                                   plane_face = "xy"
                                   plane_face = "xz"
                                   plane_face = "yz"
                                   plane_face = "xyz":
                                        where "xyz" makes level sets on all three faces.
                                   you can combine two or these x[i]x[j] (i !=j) by separatung them
                                   with "_". Example: plane_face = "xy_yz" makes contour line on the xy-plane and the yz-plane

                                   you can even switch the letters to plane_face = "zy_yx" and still get teh same result
                                   (or maybe not, i can't bother, i must find a more efficient way to code this)

                                   plane_face = None:
                                        does not make any contour lines
                          appearence = "wireframe":
                              creates a wireframe instead of a surface of the function
                          fucntion_list (list):
                               is a list of fucntions to be plotted in the same window
                    RETURNS:
                         a 3D plot of the surface
               For 2D plotting:
                    plot_style = "2D":
                         makes a plot of ONLY the level sets in the xy-plane
                    coloroption (string):
                         this is just the coloring of the level sets
                    RETURNS:
                         a 2D plot showing the level sets.
     """
     if not fx and not fx:
          def f_x(x,y, h=1.e-20):
               fx = (f(x+h, y) -f(x,y))/h
               return fx
          def f_y(x,y, h=1.e-20):
               fy = (f(x, y+h) -f(x,y))/h
               return fy
          fx = f_x
          fy = f_y


     x = np.linspace(x1,x2,N)
     y = np.linspace(y1,y2,N)
     X, Y = np.meshgrid(x,y)
     L = np.array([[np.array([X[i,j], Y[i,j]]) for i in range(X.shape[0])] for j in range(Y.shape[0])])

     def lvl_check(f, v, lvl_vals=lvl_vals,tol=tol):
          if not lvl_vals:
               pass
          else:
               for c in lvl_vals:
                   if c+tol > f(*v) > c-tol:
                        return c
               return np.inf

     if plot_style == "2D":
          if not lvl_vals:
               raise Exception("you must give some values of constant z")
          A = np.array([[lvl_check(f, L[i,j]) for i in range(N)] for j in range(N)])
          fig, ax = plt.subplots()
          ax.pcolor(X,Y,A, cmap=f"{coloroption}")

     elif plot_style == "3D":
          fig = plt.figure()
          ax = fig.gca(projection="3d")
          A = np.array([[f(*L[i,j]) for i in range(N)] for j in range(N)])
          ax.plot_surface(X,Y,A,alpha=opacity)
          if appearence == "wireframe":
               ax.plot_wireframe(X,Y,A,alpha=opacity, rstride=wireframedensity, cstride=wireframedensity)

          if not function_list and function_list != g:
               B = np.array([[g(*L[i,j]) for i in range(N)] for j in range(N)])
               if appearence == "surface":
                    ax.plot_surface(X,Y, B,alpha=opacity)

#                    slda_ax = plt.axes([0.05, 0.5, 0.3, 0.03])
#                    # axes for a
#                    sldd_ax =plt.axes([0.05, 0.4, 0.3, 0.03]) # axes for d
#                    ax = plt.axes([0.5, 0.1, 0.4, 0.8])        # axes for plot
#                    u = np.linspace(x1,x2, N)
#                    v =  np.linspace(y1,y2, N)
#                    a0, d0 = 1., 0.05
#                    slda = mpw.Slider(slda_ax, 'a', x1, x2, valinit=a0, valfmt="%1.2f")
#                    sldd = mpw.Slider(sldd_ax, 'd', y1, y2, valinit=d0,  valfmt="%1.3f")
#
#                    def L(x,y, a,b):
#                         Tplane = f(a,b) + fx(a,b)*(x-a) + fy(a,b)*(y-b)
#                         return Tplane
#
#                    P = np.array([[L(*L[i,j], u,v) for i in range(N)] for j in range(N)])
#
#                    plane = ax.plot_surface(X,Y,P)
#
#                    def update_curve1(val):
#                         a = slda.val
#                         d = sldd.val
#                         line.set_xdata(np.cos(a*u)*np.exp(-a*u))
#                         line.set_ydata(np.sin(d*v)*np.exp(-d*v))
#                    #     line.set_xdata(p1(a,d) + p1x(a,d)*(u-a) p1y(a,d)*(v-d))
#
#                         slda.on_changed(update_curve1)
#                         sldd.on_changed(update_curve1)



          elif function_list != None and len(function_list) != 0:
               L = [np.array([[function(*L[i,j]) for i in range(N)] for j in range(N)]) for function in function_list]
               for item in L:
                    ax.plot_surface(X,Y,item,alpha=opacity)


          if not plane_face:
               pass
          else:
               #plot contour projection of the corespondingly fed axis face
               if plane_face == "xy":
                    ax.contour(X,Y, A, zdir="z", offset = -(np.abs(min(lvl_vals))+ np.abs(max(lvl_vals)))/2)
               elif plane_face == "xz":
                    ax.contour(X,Y, A, zdir="y", offset = (np.abs(y1) + np.abs(y2))/2)
               elif plane_face == "yz":
                    ax.contour(X,Y, A, zdir="x", offset = -(np.abs(x1) + np.abs(x2))/2)
               elif plane_face == "xy_xz" or plane_face =="xz_xy":
                    ax.contour(X,Y, A, zdir="z", offset = -(np.abs(min(lvl_vals))+ np.abs(max(lvl_vals)))/2)
                    ax.contour(X,Y, A, zdir="y", offset = (np.abs(y1) + np.abs(y2))/2)
               elif plane_face == "xy_yz" or plane_face == "yz_xy":
                    ax.contour(X,Y, A, zdir="z", offset = -(np.abs(min(lvl_vals))+ np.abs(max(lvl_vals)))/2)
                    ax.contour(X,Y, A, zdir="x", offset = -(np.abs(x1) + np.abs(x2))/2)
               elif plane_face == "yz_xz" or plane_face == "xz_yz":
                    ax.contour(X,Y, A, zdir="x", offset = -(np.abs(x1) + np.abs(x2))/2)
                    ax.contour(X,Y, A, zdir="y", offset = (np.abs(y1) + np.abs(y2))/2)
               elif plane_face == "xyz":
                    ax.contour(X,Y, A, zdir="x", offset = -(np.abs(x1) + np.abs(x2))/2)
                    ax.contour(X,Y, A, zdir="z", offset = -(np.abs(min(lvl_vals))+ np.abs(max(lvl_vals)))/2)
                    ax.contour(X,Y, A, zdir="y", offset = (np.abs(y1) + np.abs(y2))/2)
               else:
                    raise Exception("you have to choose a value for the paramater plane_face, refer to help for more info")
          #set limits
          ax.set_xlim3d(x1,x2)
          ax.set_ylim3d(y1,y2)
          ax.set_zlim3d(np.min(A), np.max(A))
          #set labels
          ax.set_xlabel("X axis")
          ax.set_ylabel("Y axis")
          ax.set_zlabel("Z axis")



fig = plt.figure()
slda_ax = plt.axes([0.05, 0.5, 0.3, 0.03])
# axes for a
sldd_ax =plt.axes([0.05, 0.4, 0.3, 0.03]) # axes for d
ax = plt.axes([0.5, 0.1, 0.4, 0.8])        # axes for plot
u = np.linspace(0, 14*np.pi, 400)
v =  np.linspace(0, 14*np.pi, 400)
a0, d0 = 1., 0.05
slda = mpw.Slider(slda_ax, 'a', 0, 2, valinit=a0, valfmt="%1.2f")
sldd = mpw.Slider(sldd_ax, 'd', 0, 0.1, valinit=d0,  valfmt="%1.3f")

ax.set_xlim([-1, 1])
ax.set_ylim([-1, 1])
#ax.set_zlim([-1, 1])
ax.set_aspect(1.0)
# the aspect ratio should be 1 here
line = ax.plot(np.cos(a0*t)*np.exp(-d0*t), np.sin(t)*np.exp(-d0*t))
def update_curve1(val):
     a = slda.val
     d = sldd.val
     line.set_xdata(np.cos(a*u)*np.exp(-a*u))
     line.set_ydata(np.sin(d*v)*np.exp(-d*v))
#     line.set_xdata(p1(a,d) + p1x(a,d)*(u-a) p1y(a,d)*(v-d))

     slda.on_changed(update_curve1)
     sldd.on_changed(update_curve1)



Q = [p2, p3]
A = list(np.linspace(-1,1,10))
level_sets(p1, -20,20,-20,20,200, plot_style="3D", function_list=Q)
























def level_sets1(f,lvl_vals, x1,x2, y1,y2, N, tol=1.e-2):
     X=list(np.linspace(x1,x2, N))
     Y=list(np.linspace(y1,y2, N))
     L = [([],[]) for c in lvl_vals]
     for ind, c in enumerate(lvl_vals):
          for i in X:
               for j in Y:
                    new = f(i,j)
                    if c+tol > new > c-tol:
                         L[ind][0].append(i)
                         L[ind][1].append(j)
     plt.figure()
     for ind,i in enumerate(L):
          plt.scatter(i[0], i[1], label = f"f(x,y) = {lvl_vals[ind]}")
     plt.show()

level_sets(p, range(-5,5), -10,10,-10,10,1000)





def f(x):
     return x**2
def g(x):
     return -(x-1)/(x+1)




X = np.linspace(-1000,-1, 1000)
Y = np.linspace(1,1000,1000)
a = [f(X) for i in X]
b = [g(X) for i in X]
plt.figure()
plt.plot(X, g(X), label = "transformed lower bounds")
plt.figure()
plt.plot(Y, g(Y), label = "transformed lower bounds")
plt.show()









print(np.arccos(3/5), np.arcsin(4/5))


def direct_para_normal(p=None,u=None,v=None, point1=None, point2=None, point3=None):

     def motor(p,u,v):
          p,u,v = np.array(p), np.array(u), np.array(v)
          N = np.cross(v,u)
          N[1] = -N[1]
          n = np.array([-(u[1]*N[1] + u[2]*N[2]), u[0]*N[1], u[0]*N[2]])
          a = np.array([np.matmul(n,p)])
          return np.hstack([n,a])

     if not point1 and point1 == point2==None:
          if not p and p ==u==v:
               raise Exception("a plane has to be define with at least three points or a point and two non-colinear vectors")
          return motor(p,u,v)
     else:
          p1, p2, p3 = point1, point2, point3
          u = p2 - p1
          v = p3 - p1
          return motor(p1, u,v)
#print(direct_para_normal(p=[0,0,0,0.], u=[-2,1,0,1.], v=[5,-4,1,0.]))
#print(direct_para_normal(point1=[1,1,0.], point2=[0,0,0.], point3=[2,0,1.]))







def many_iter(M):
     """
     for this function you havr to manually apply the function each time to
     your matrix, and then apply slcicing M[1:M.shape[0]-i, 1:M.shape[1]-i]
     where i depends on how many iterations you have done. i starts with a value of zero.
     example:
          T = np.array([[1,6,3,0.],
                       [1,9,1,1],
                       [5,5,7,5]])

          T1 = many_iter(T)
          T2 = T1[1:3, 1:4]

          T3 = many_iter(T2)
          print(T)


          >>> [[  1.   6.   3.   0.]
               [  0.   3.  -2.   1.]
               [  0.   0. -74.  40.]]

     in this exaplme the function many_iter is applied twice
     and we're done. the result is a triangular matrix. with which
     you now can find a basis for the null-space using simple subtitution.

     """
     def subtraction(M):

          a,b  = M.shape[0], M.shape[1]

          M1 = M[1:a, 0:b]*M[0,0]
          M2 = np.array([M[0]*M[i,0] for i in range(1, a)])
          M[1:a, 0:b]  = M1 - M2
          return M

     lit = [M]
     for i in range(1):
          new = subtraction(lit[-1])
          lit.append(new)
     return lit[-1]

K = np.array([[1,1,1.],
              [2,-1,5],
              [3,2,4],
              [4,3,5]])
K1 = many_iter(K)
K2 = K1[1:4, 1:3]
K3 = many_iter(K2)
print(K)

K = np.array([[1,1,2,1,2.],
             [0,1,1,-1,2],
             [2,1,-1,2,-1],
             [1,-2,-2,1,1]])
K1 = many_iter(K)
K2 = K1[1:4, 1:5]
K3 = many_iter(K2)
K4 = K2[1:3, 1:4]
K5 = many_iter(K4)
K6 = K5[1:2, 1:3]
print(K)












K = np.array([[1,1,2,5.],
              [3,2,1,1],
              [2,3,2,3],
              [1,4,3,5]])
K1 = many_iter(K)
K2 = K1[1:4, 1:4]
K3 = many_iter(K2)
K4 = K3[1:3, 1:3]
print(K)







def f(M):
     a,b  = M.shape[0], M.shape[1]
     def subtraction(n,M):
          def function(M):
               """
               Ö should always start with 0
               """
               M1 = M[1:a, 0:b]*M[0,0]
               M2 = np.array([M[0]*M[i,0] for i in range(1, a)])
               M[1:a, 0:b]  = M1 - M2
               return M

          lit = [M]
          for i in range(1):
               new = function(lit[-1])
               lit.append(new)
          p,q = lit[-1].shape[0], lit[-1].shape[1]
          second = lit[-1][1:p-n, 1:q-n]
          return second

     list_ = [M]
     for k in range(0, b):
          new = subtraction(k,list_[-1])
          list_.append(new)
          print(M)
     return M




K = np.array([[1,3,4,1.],
              [1,1,2,-1],
              [2,5,7,2]])
np.max(K)
print(many_iter(K))


K1 = many_iter(0,K)
K3 = many_iter(1,K1)
print(K)




E = np.array([[5,1,2.,1],
              [1,2,1, 1],
              [3,3,2, 1],
              [5,4,3, 2]])


E1 = many_iter(E)
E2 = E1[1:4, 1:3]
E3 = many_iter(E2)
E4 = E3[1:3, 1:2]
E5 = many_iter(E4)
E6 = E5[1:2, 1:1]

print(E)




A = np.array([[1,2,1.],
              [2,3,1],
              [3,2,1],
              [1,3,2],
              [2,1,3]])

e = np.linalg.eigvals()
A1 = many_iter(A)
A2 = A1[1:5, 1:3]
A3 = many_iter(A2)

print(A)




A = np.array([[1,2,1.],
             [2,3,1],
             [3,2,-1],
             [1,3,2]])
A1 = many_iter(A)
A2 = A1[1:4, 1:3]
A3 = many_iter(A2)

print(A)















P = np.array([[1,2,3,4,5.,],
              [1,3,4,5,6],
              [2,5,7,9,11]])
P1 = many_iterq(P)
print(P)

K = np.array([[3,2,1.,39],
              [2,3,1, 34],
              [1,2,3,26]])
K1 = many_iter(K)
K2 = K1[1:3, 1:4]
k3 = many_iter(K2)
print(K)




T = np.array([[1,6,3,0.],
             [1,9,1,1],
             [5,5,7,5]])
T1 = many_iter(T)
T2 = T1[1:3, 1:4]

T3 = many_iter(T2)
print(T)




A = np.array([[1,2,-1,0,4,5,7,8],
              [9,1,0,-1.,6,3,6,7],
              [7,5,3,4,9,6,2,6,],
              [8,4,7,3,6,8,3,8],
              [10,6,8,3,7,9,3,78]])

E1 = many_iter(A)
E2 = E1[1:5, 1:8]
E3 = many_iter(E2)
E4 = E3[1:4, 1:7]
E5 = many_iter(E4)
E6 = E5[1:3, 1:6]
print(E6)


q = E[1:3, 1:4]
s = many_iter(0,0,q)
print(s)
print(E)
print(q)
print(E)
print(many_iter(1,1, E))



A1 = np.array([A[1]*A[i,1] for i in range(2, A.shape[0])])
A2 = np.row_stack([np.zeros_like(A[1]), A1])
print(A2)
print(A[2:4, 1:3])



R = np.array([1,2,3,4,5,6,7,8,9])
print(R[1:3])



print(np.linalg.inv(A))



def lin_dep(Set, margin=1.e-5, error_notification="raise_excpetion"):
    """
    This function checks if a set of vectors are linearly dependent or not.

    One of the conditions that makes a set of vectors linearly independent is that Xa = 0
    results in a vector a whose components all are zero.
    if they would have been anything else then they would be linearly dependent

    On INPUTS:
        Set (list):
            a list of arrays
        error_notification (string/bool):
            can either be the string "raise_excpetion" or the boolena  value None
    ON RETURN:


    """
    Dim = []
    for i, s in enumerate(Set):
        if type(s) != np.ndarray or np.ndim(s) != 1:
            if not error_notification:
                Set = Set[:i] + Set[-(len(Set) -(i+1)):]
            else:
                raise Exception("Set must only contain vectors")

    for s in Set:
        Dim.append(s.shape[0])

    if max(Dim) == min(Dim) < len(Dim):
        return "linearly dependent 0: amount of vectors is larger than their dimension"
    else:
        m_Dim = max(Dim)
        ind = Dim.index(m_Dim)
        M = np.zeros((len(Set), Set[ind].shape[0]))
        for i, s in enumerate(Set):
            s_z = np.zeros(m_Dim - Dim[i])
            s = np.hstack([s, s_z])
            M[i] = s

        zero = np.zeros(m_Dim)
        try:
            conclusion = np.linalg.solve(M, zero)
        except sl.LinAlgError:
            return np.linalg.det(M),"linearly dependent 1", True

        if np.linalg.norm(conclusion) < margin:
            return f"det = {np.linalg.det(M)}",conclusion, "linearly independent 2", False
        else:
            return M,conclusion, "linearly dependent 3", True


K = np.array([[1,0,0,0],
              [1,-1,0,2],
              [1,0,-1,0],
              [0,2,-1,1]])
Set = list(K)
print(lin_dep(Set))




print(np.array([0,0,0,0]) + np.array([1,1]))

#E = np.array([[0., 1.,  2.,  33.],
#              [1., 5.,  6.,  66.],
#              [2., 9.,  10., 84.],
#              [3., 12., 13., 99.],
#              [4., 15., 16., 46.],
#              [5., 18., 19., 87.],
#              [6., 21., 22., 45.],
#              [7., 24., 25., 77.],
#              [8., 27., 28., 44.]])
#print(E[1:,])






def angle_test(M, margin=1.e-9):
     theta = np.array([[np.arccos(np.matmul(M[i],M[j])/(np.linalg.norm(M[i])*np.linalg.norm(M[j]))) for i in range(M.shape[0])] for j in range(M.shape[0])])

     return theta


print(angle_test(A.T))




def machine1(n,r):
     for j in range(r.shape[1]-1):
          for i in range(n,r.shape[1]-1):
               r[i] = r[j-1,j-1]*r[i]
               r[i] = r[i] - r[j,j-1]*r[i-1]
     return r

print(A)
print(machine1(1,A))





































def plane_intersect_plane(n,a, m,b, M1_para=None, M2_para=None):
     """
     this is an upgraded version of the plane intersection function. this one
     can even handle planes given in parametric form.
     the parametric form is given with a matrix Mi_para containing
     a pointin 3D space and two vectors in that order.



     """
     def motor1(n,a, m,b):
          if type(n) != np.ndarray != type(m):
               if type(n) != list != type(m)  or type(n) != tuple != type(m):
                    raise Exception("n and m must be vectors or lists,tuple that can be turned to vectors")
               else:
                    n=np.array(n)
                    m=np.array(m)
                    if n.shape[0] != m.shape[0]:
                         raise Exception("n and m must be of the same dimension")
          q = b*n[0] - a*m[0]
          M = np.array([[0, -n[2], n[1]],[n[2], 0, -n[0]],[-n[1], n[0], 0]])
          N = np.matmul(M,m)
          N[1] = -N[1]
          line_direction_vector = np.array([n[1]/n[0] - (n[2]*N[2])/(n[0]*N[1]), -1, N[2]/N[1]])
          line_point = np.array([(a - q/N[2])/n[0], q/N[2], 0])
          return np.array([line_point,line_direction_vector])

     if not M1_para and M1_para == M2_para:
          return motor1(n,a, m,b)

     else:
          if type(M1_para) != np.ndarray != type(M2_para):
               if type(M1_para) != list != type(M2_para)  or type(M1_para) != tuple != type(M2_para):
                    raise Exception("M1_para and M2_para must be matrices or lists,tuple that can be turned to matices")
               else:
                    M1_para=np.array(M1_para)
                    M2_para=np.array(M2_para)
                    if M1_para.shape != M1_para.shape:
                         raise Exception("M1_para and M2_para must be of the same dimension")
          def para_normal(p,u,v):
               """
               thsi function take a point (p) and two vectors (u,v) that are used for
               the equation for a plane in parametric form to produce the normal vector
               of that plane
               """
               N = np.cross(v,u)
               N[1] = -N[1]
               e = np.array([-(u[1]*N[1] + u[2]*N[2]), u[0]*N[1], u[0]*N[2]])
               return e
          if not m and m == b:
               if not M1_para and M1_para != M2_para:
                    m = para_normal(*M2_para)
                    b = np.matmul(M2_para[0], m)
               elif not M2_para and M1_para != M2_para:
                    m = para_normal(*M2_para)
                    b = b = np.matmul(M2_para[0], m)

          if not n and n ==m==a==b:
               m = para_normal(*M2_para)
               b = np.matmul(M2_para[0], m)
               n = para_normal(*M1_para)
               a = b = np.matmul(M1_para[0], n)

          return motor1(n,a,m,b)

print(plane_intersect_plane(np.array([2,1,-3.]),7,np.array([1,2.,-1]),4))




E = np.array([[0., 1.,  2.,  33.],
              [1., 5.,  6.,  66.],
              [2., 9.,  10., 84.],
              [3., 12., 13., 99.],
              [4., 15., 16., 46.],
              [5., 18., 19., 87.],
              [6., 21., 22., 45.],
              [7., 24., 25., 77.],
              [8., 27., 28., 44.]])
print(*E)


def f(z, b=None, c=None):
     if not b and b==c :
          return 1
     else:
          return z
print(f(3))





















def shortest_distance_line_extPoint(external_point, line_point1=None, line_point2=None, line_start_point=None, line_direction_vector=None, margin=1.e-10):
     """
     This fucntion calculates the shortest distance between a line and an externeal point

     INputs:
          external_point (array)

          optional inputs:
               line_point1, line_point2 (array):
                    are there to define a line in space to work with
               or

               line_start_point (array):
                    this poibtis given in any parameterised representation of a straight line in space
                    line_direction_vector (array):
                    this is the vector giving the direction you must go to find points on the line

          Observe that inly one of the two sets of optional inputs must have a value of None.
     OUTput:
          LEN_QR (float):
               is the resulting length between the line and the external point
     """
     A,B,R = line_point1, line_point2, external_point
     if not line_direction_vector and not line_start_point:
          if not A and not B:
               raise Exception("you forgot to defineline_point1 and line_point2")

          if type(A) != np.ndarray != type(B) != type(R):
               if type(A) != list != type(B) != type(R) or type(A) != tuple != type(B) != type(R):
                    raise Exception("point_line1, and point_line2 must be vectors or lists,tuple that can be turned to vectors")
          else:
               A=np.array(A)
               B=np.array(B)
               R=np.array(R)

               if A.shape[0] != B.shape[0] != R.shape[0]:
                    raise Exception("line_point1, line_point2, external_point must be of same dimension")
          u = B - A

     elif not line_point1 and not line_point2:
          A = line_start_point
          u = line_direction_vector
          if type(A) != np.ndarray != type(u) != type(R):
               if type(A) != list != type(u) != type(R) or type(A) != tuple != type(u) != type(R):
                    raise Exception("point_line1, and point_line2 must be vectors or lists,tuple that can be turned to vectors")
          else:
               A=np.array(A)
               u=np.array(u)
               R=np.array(R)

               if A.shape[0] != u.shape[0] != R.shape[0]:
                    raise Exception("line_point1, line_point2, external_point must be of same dimension")


     AR = R - A
     T = np.matmul(AR, u)/np.linalg.norm(u)**2
     QR = AR - u*T
     Len_QR = np.linalg.norm(QR)
     return Len_QR

print(shortest_distance_line_extPoint([1.,2,1],line_start_point=[-2.,1,1],line_direction_vector=[2.,-2,1]))




def plane_intersect_plane(n,a, m,b):
     """
     this fucntion calculates
     """
     if type(n) != np.ndarray != type(m):
          if type(n) != list != type(m)  or type(n) != tuple != type(m):
               raise Exception("point_line1, and point_line2 must be vectors or lists,tuple that can be turned to vectors")
          else:
               n=np.array(n)
               m=np.array(m)
               if n.shape[0] != m.shape[0]:
                    raise Exception("line_point1, line_point2, external_point must be of same dimension")
     q = b*n[0] - a*m[0]
     M = np.array([[0, -n[2], n[1]],[n[2], 0, -n[0]],[-n[1], n[0], 0]])
     N = np.matmul(M, m)
     line_direction_vector = np.array([n[1]/n[0] - (n[2]*N[2])/(n[0]*N[1]), -1, N[2]/N[1]])
     line_point = np.array([(a - q/N[2])/n[0], q/N[2], 0])

     return np.array([line_direction_vector, line_point])

print(plane_intersect_plane(np.array([5,3,1.]),3,np.array([1,1,1.]),0))



def many_iter(M):
     """
     for this function you havr to manually apply the function each time to
     your matrix, and then apply slcicing M[1:M.shape[0]-i, 1:M.shape[1]-i]
     where i depends on how many iterations you have done. i starts with a value of zero.
     example:
          T = np.array([[1,6,3,0.],
                       [1,9,1,1],
                       [5,5,7,5]])

          T1 = many_iter(T)
          T2 = T1[1:3, 1:4]

          T3 = many_iter(T2)
          print(T)

          >>> [[  1.   6.   3.   0.]
               [  0.   3.  -2.   1.]
               [  0.   0. -74.  40.]]

     in this exaplme the function many_iter is applied twice
     and we're done. the result is a triangular matrix. with which
     you now can find a basis for the null-space using simple subtitution.

     """
     def subtraction(M):
          """
          Ö should always start with 0
          """
          a,b  = M.shape[0], M.shape[1]

          M1 = M[1:a, 0:b]*M[0,0]
          M2 = np.array([M[0]*M[i,0] for i in range(1, a)])
          M[1:a, 0:b]  = M1 - M2
          return M

     lit = [M]
     for i in range(1):
          new = subtraction(lit[-1])
          lit.append(new)
     return lit[-1]




def lin_dep(Set, margin=1.e-5, error_notification="raise_excpetion"):
    """
    This function checks if a set of vectors are linearly dependent or not.

    One of the conditions that makes a set of vectors linearly independent is that Xa = 0
    results in a vector a whose components all are zero.
    if they would have been anything else then they would be linearly dependent

    On INPUTS:
        Set (list):
            a list of arrays
        error_notification (string/bool):
            can either be the string "raise_excpetion" or the boolena  value None
    ON RETURN:


    """
    Dim = []
    for i, s in enumerate(Set):
        if type(s) != np.ndarray or np.ndim(s) != 1:
            if not error_notification:
                Set = Set[:i] + Set[-(len(Set) -(i+1)):]
            else:
                raise Exception("Set must only contain vectors")

    for s in Set:
        Dim.append(s.shape[0])

    if max(Dim) == min(Dim) < len(Dim):
        return "linearly dependent 0: amount of vectors is larger than their dimension"
    else:
        m_Dim = max(Dim)
        ind = Dim.index(m_Dim)
        M = np.zeros((len(Set), Set[ind].shape[0]))
        for i, s in enumerate(Set):
            s_z = np.zeros(m_Dim - Dim[i])
            s = np.hstack([s, s_z])
            M[i] = s

        zero = np.zeros(m_Dim)
        try:
            conclusion = np.linalg.solve(M, zero)
        except sl.LinAlgError:
            return np.linalg.det(M),"linearly dependent 1", True

        if np.linalg.norm(conclusion) < margin:
            return f"det = {np.linalg.det(M)}",conclusion, "linearly independent 2", False
        else:
            return M,conclusion, "linearly dependent 3", True

S = [np.array([1,0,1,0.]),np.array([1,1,1,1.]), np.array([1,2,2,1.]), np.array([3,4,1,3])]
print(lin_dep(S))























#x, y = sp.symbols('x y', real=True)
#f1 = x**8-28*x**6*y**2+70*x**4*y**4+15*x**4-28*x**2*y**6-90*x**2*y**2+y**8+15*y**4-16
#f2 = 8*x**7*y-56*x**5*y**3+56*x**3*y**5+60*x**3*y-8*x*y**7-60*x*y**3
#print((diff(f1, x),diff(f1,y)))
#print((diff(f2,x), diff(f2,y)))



def Complete(F,J, a,b,c,d,N, coloroption, margin=1.e-4, plot="iter"):
     zero =[]
     def Newt_motor(F,J, v0,max_iter=100, margin=margin, step=0.0001, h=1.e-9, finite=False, once=False):
          xvals   = [v0]
          x_max   = 1/margin
          k       = 0

          def solving_motor(v):
              def univ_finite_diff(v, h):
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
                      delta = np.linalg.solve(univ_finite_diff(v, h), -F(v))
                      new = delta + v
                      return new
                  elif once == True:
                      J_once = univ_finite_diff(v0, h)
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

     def zero_collect_motor(F, J, v):
          """
          Task 2
          """
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

     x = np.linspace(a,b,N)
     y = np.linspace(b,c,N)
     X, Y = np.meshgrid(x,y)
     L = np.array([[np.array([X[i,j], Y[i,j]]) for i in range(X.shape[0])] for j in range(Y.shape[0])])
     if plot == "zero":
          A = np.array([[zero_collect_motor(F, J, L[i,j])[0] for i in range(N)] for j in range(N)])
     elif plot == "iter":
           A = np.array([[zero_collect_motor(F, J, L[i,j])[1] for i in range(N)] for j in range(N)])

     fig, ax = plt.subplots()
     ax.pcolor(X,Y,A, cmap=f"{coloroption}")








def Test_F(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([3*x**3-3*x*y**2 -1 , 3*x**2*y -y**3])



def Test_J(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[(9*x**2 - 3*y**2), -6*x*y],
                     [6*x*y, 3*x**2 - 3*y**2]])


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])



def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])

def Test_F2(X):
    X = np.array(X)
    return np.array([X[0]**3 - 3*X[0]*X[1]**2 - 2*X[0] -2, 3*X[0]**2*X[1] - X[1]**2 -2*X[1]])


def Test_J2(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2 -2 , -6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 2*X[1] -2]])


def Test_F3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([x**8-28*x**6*y**2+70*x**4*y**4+15*x**4-28*x**2*y**6-90*x**2*y**2+y**8+15*y**4-16,
                     8*x**7*y-56*x**5*y**3+56*x**3*y**5+60*x**3*y-8*x*y**7-60*x*y**3])


def Test_J3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]

    return np.array([[8*x**7 - 168*x**5*y**2 + 280*x**3*y**4 + 60*x**3 - 56*x*y**6 - 180*x*y**2, -56*x**6*y + 280*x**4*y**3 - 168*x**2*y**5 - 180*x**2*y + 8*y**7 + 60*y**3 ],
                     [56*x**6*y - 280*x**4*y**3 + 168*x**2*y**5 + 180*x**2*y - 8*y**7 - 60*y**3, 8*x**7 - 168*x**5*y**2 + 280*x**3*y**4 + 60*x**3 - 56*x*y**6 - 180*x*y**2]])




Complete(Test_F2, Test_J2,-2,2,-2,2,500, "plasma")
#Complete(Test_F2, Test_J2,-2,2,-2,2,500, "PuBuGn")

plt.show()


























































def motor_integration(F, J, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, coloroption="plasma"):
     """
     This Fucntion conatains all the other Tasks and methods  integrated into one function.
     """

     zero = []
     N = int(N)

     def Newt_motor(F, J, v0, max_iter=max_iter, margin=1.e-4, step=step, h=1.e-9, finite=False, once=False):
          """
          Task 2
          """
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






     def zero_collect_motor(F, J, v):
          """
          Task 2
          """
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






#            elif A_method == "abs":
#                A = np.zeros((N,N))
#                L = tensor_motor(a,b,c,d,N)
#                for i in range(N):
#                    for j in range(N):
#                        new = Newt_motor(L[j,i])[0]
#                        A[i,j] = np.abs(new[0]) + np.abs(new[1])


#     if PLOT == True:
#          x = np.linspace(a,b,N)
#          y = np.linspace(c,d,N)
#          X,Y = np.meshgrid(x,y)
#          fig, ax = plt.subplots()
#          ax.pcolor(X, Y, A, cmap='plasma')


     def tensor_motor(a,b,c,d,N, x, y):
          """
          Part of Task 4
          """
          x = np.linspace(a,b,N)
          y = np.linspace(c,d,N)
          X,Y = np.meshgrid(x,y)
          L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
          return  L, X, Y
     L = tensor_motor(a,b,c,d,N)
     A = np.array([[zero_collect_motor(L[i,j])[1] for i in range(N)] for j in range(N)])



     X,Y = np.meshgrid(x,y)
     fig, ax = plt.subplots()
     ax.pcolor(X, Y,A, cmap=f"{coloroption}")

     return zero













def Test_F(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([3*x**3-3*x*y**2 -1, 3*x**2*y -y**3])



def Test_J(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[9*x**2 - 3*y**2, -6*x*y],
                     [6*x*y, 3*x**2 - 3*y**2]])


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])



def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])

def Test_F2(X):
    X = np.array(X)
    return np.array([X[0]**3 - 3*X[0]*X[1]**2 - 2*X[0] -2, 3*X[0]**2*X[1] - X[1]**2 -2*X[1]])


def Test_J2(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2 -2 , -6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 2*X[1] -2]])


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




print(motor_integration(-2,2,-2,2,20, Test_F1, Test_J1,coloroption="plasma"))










