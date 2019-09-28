# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 13:15:27 2019

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

def motor_integration(F, J, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT="both", A_method="lol", coloroption="plasma"):
   """
   This Fucntion conatains all the other Tasks and methods  integrated into one function.
   """

   zero = []

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

   def Matrix_A(A_method=A_method):
       """
       Task 4
       """
       def tensor_motor(a,b,c,d,N):
           """
           Part of Task 4
           """
           x = np.linspace(a,b,N)
           y = np.linspace(c,d,N)
           X,Y = np.meshgrid(x,y)
           L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
           return  L

       if A_method == "lol":
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

       A = np.zeros((int(N),int(N)))
       L = tensor_motor(a,b,c,d,N)
       A = [[zero_collect_motor(L[i,j])[0] for i in range(N)] for j in range(N)]



#            elif A_method == "abs":
#                A = np.zeros((N,N))
#                L = tensor_motor(a,b,c,d,N)
#                for i in range(N):
#                    for j in range(N):
#                        new = Newt_motor(L[j,i])[0]
#                        A[i,j] = np.abs(new[0]) + np.abs(new[1])
       return A

   if PLOT == True:
       x = np.linspace(a,b,N)
       y = np.linspace(c,d,N)
       X,Y = np.meshgrid(x,y)
       fig, ax = plt.subplots()
       ax.pcolor(X, Y, Matrix_A(A_method), cmap='plasma')

   elif PLOT == False:
       return Matrix_A(A_method)

   elif PLOT =="both":
       x = np.linspace(a,b,N)
       y = np.linspace(c,d,N)
       X,Y = np.meshgrid(x,y)
       fig, ax = plt.subplots()
       ax.pcolor(X, Y,Matrix_A(A_method), cmap=f"{coloroption}")

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



#print(motor_integration(-2,2,-2, 2, 10, Test_F2, Test_J2, PLOT="both",plotting_method="mine", coloroption="Spectral"))
#print(motor_integration(-2,2,-2, 2,10, Test_F2, Test_J2, PLOT="both",plotting_method="fml", coloroption="gray"))
print(motor_integration(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both", coloroption="plasma"))
#print(motor_integration(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both", coloroption="gray"))
#print(motor_integration(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both", coloroption="plasma"))
#print(motor_integration(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both", coloroption="gray"))
#print(motor_integration(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both", coloroption="plasma"))
#print(motor_integration(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both", coloroption="gray"))






























































def Unified_PLOT(F, J, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT="both", A_method="int", coloroption="plasma"):
        """
        This Fucntion conatains all the other Tasks and methods  integrated into one function.
        """

        zero = []

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

        def Matrix_A(A_method):
            """
            Task 4
            """
            def tensor_motor(a,b,c,d,N):
                """
                Part of Task 4
                """
                x = np.linspace(a,b,N)
                y = np.linspace(c,d,N)
                X,Y = np.meshgrid(x,y)
                L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
                return  L

            #if A_method == "int":
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

        A = np.zeros((N,N))
        L = tensor_motor(a,b,c,d,N)
        A = [[zero_collect_motor(L[i,j])[1] for i in range(N)] for j in range(N)]



#            elif A_method == "abs":
#                A = np.zeros((N,N))
#                L = tensor_motor(a,b,c,d,N)
#                for i in range(N):
#                    for j in range(N):
#                        new = Newt_motor(L[j,i])[0]
#                        A[i,j] = np.abs(new[0]) + np.abs(new[1])


        if PLOT == True:
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, A, cmap='plasma')

        elif PLOT == False:
            return Matrix_A(A_method)

        elif PLOT =="both":
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
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



#print(motor_integration(-2,2,-2, 2, 10, Test_F2, Test_J2, PLOT="both",plotting_method="mine", coloroption="Spectral"))
#print(motor_integration(-2,2,-2, 2,10, Test_F2, Test_J2, PLOT="both",plotting_method="fml", coloroption="gray"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both",A_method="int", coloroption="plasma"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both",A_method="int", coloroption="gray"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both",A_method="int", coloroption="plasma"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both",A_method="int", coloroption="gray"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both",A_method="int", coloroption="plasma"))
print(Unified_PLOT(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both",A_method="int", coloroption="gray"))



def Unified_PLOT(self, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT="both", A_method="int", coloroption="plasma"):
        """
        This Fucntion conatains all the other Tasks and methods  integrated into one function.
        """
        F = self.F
        J = self.J
        zero = self.zero

        def Newt_motor(self, v0, max_iter=max_iter, margin=1.e-4, step=step, h=1.e-9, finite=False, once=False):
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
                    self.J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                                       [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                    return self.J

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

        def Matrix_A(A_method):
            """
            Task 4
            """
            def tensor_motor(a,b,c,d,N):
                """
                Part of Task 4
                """
                x = np.linspace(a,b,N)
                y = np.linspace(c,d,N)
                X,Y = np.meshgrid(x,y)
                L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
                return  L

            if A_method == "int":
                def zero_collect_motor(F, J, v):
                    """
                    Task 2
                    """
                    v_a, v_i = Newt_motor(self, v)
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

                A = np.zeros((N,N))
                L = tensor_motor(a,b,c,d,N)
                A = [[zero_collect_motor(L[i,j])[0] for i in range(N)] for j in range(N)]



            elif A_method == "abs":
                A = np.zeros((N,N))
                L = tensor_motor(a,b,c,d,N)
                for i in range(N):
                    for j in range(N):
                        new = Newt_motor(L[j,i])[0]
                        A[i,j] = np.abs(new[0]) + np.abs(new[1])
            return A

        if PLOT == True:
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, Matrix_A(A_method), cmap='plasma')

        elif PLOT == False:
            return Matrix_A(A_method)

        elif PLOT =="both":
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y,Matrix_A(A_method), cmap=f"{coloroption}")

            return zero



























class Fractal2D_amer:
    """
    Task 1
    """
    def __init__(self, F, J):
        self.F = F
        self. J = J
        self.zero = []

    def Complete(self, a,b,c,d,N, coloroption="plasma", margin=1.e-4, plot="iter"):
         zero = self.zero
         F = self.F
         J = self.J
         def Newt_motor(self, v0, max_iter=100, margin=margin, step=0.0001, h=1.e-9, finite=False, once=False):
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




class fractal2D_erik:
    def __init__(self, f, j=None):
        self.f=f
        self.zeroes = [] #our list of zeroes
        self.j=j
        if j is None:
            self.j = lambda x, y: self.finite_difference(f,x,y)
    def finite_difference(self, f, x, y, h=1.e-9):
        return np.array([[(self.f(x+h,y)[0] -self.f(x,y)[0])/h, (self.f(x,y+h)[0] -self.f(x,y)[0])/h],
                         [(self.f(x+h,y)[1] -self.f(x,y)[1])/h, (self.f(x,y+h)[1] -self.f(x,y)[1])/h]])

    def newton(self, x0, maxiteration=1.e6, tol=1.e-4):
        x = np.array(x0, dtype=float) #makes it an
        f_value=self.f(*x) #unzips the array and initialize it with our function
        f_norm = np.linalg.norm(f_value) #calculates our norm(length)
        iteration_counter = 0
        x_max = 1/tol
        while abs(f_norm) > tol and iteration_counter < maxiteration and np.all(x<abs(x_max)):
            delta = np.linalg.solve(self.j(*x), -f_value) #solves for delta
            x += delta
            f_value = self.f(*x)
            f_norm = np.linalg.norm(f_value)
            iteration_counter += 1
        if abs(f_norm) > tol: #return -1 if the length is still not close to zero
            iteration_counter = -1
        return x,iteration_counter #returns the coordinate and iteration

    def simplified_newton(self, x0, maxiteration=1.e6, tol=1.e-4):
        x = np.array(x0, dtype=float) #makes it an array
        f_value=self.f(*x) #unzips the array and initialize it with our function
        f_norm = np.linalg.norm(f_value) #calculates our norm(length)
        jacobian = self.j(*x) #only difference is that we calculate jacobian one time with x0 and never updates it
        iteration_counter = 0
        x_max = 1/tol
        while abs(f_norm) > tol and iteration_counter < maxiteration and np.all(x<abs(x_max)):
            delta = np.linalg.solve(jacobian, -f_value) #solves for delta
            x += delta
            f_value = self.f(*x)
            f_norm = np.linalg.norm(f_value)
            iteration_counter += 1
        if abs(f_norm) > tol: #return -1 if the length is still not close to zero
            iteration_counter = -1
        return x, iteration_counter#returns the coordinate

    def zeroes_index(self, x0, simplified, tol=1.e-4):
        if not simplified:
            res,iteration_counter = self.newton(x0)
        else:
            res,iteration_counter = self.simplified_newton(x0)
        if (res==np.inf).any(): #check infinity
            return (np.inf,iteration_counter)
        if not len(self.zeroes): #check if list is empty
           self.zeroes.append(res)
           return (0,iteration_counter)
        for k,i in enumerate(self.zeroes): #k is the iteration nr, i is the list value
            if np.linalg.norm(res-i) < tol: #check the lenght is less than tol
                return (k,iteration_counter)
        self.zeroes.append(res)
        return (len(self.zeroes)-1,iteration_counter) #needs to be -1 here because we append before the return

    def plot(self, N, a, b, c, d, simplified):
        x = np.linspace(a,b,N)#a,b defines the x-axis length and N how many points inbetween
        y = np.linspace(c,d,N)
        X, Y = np.meshgrid(x,y) #stores a grid in 2 matrices
        A = np.array([np.array([self.zeroes_index(np.array([X[i][j],Y[i][j]]), simplified)[0] for i in range (N)]) for j in range(N)])
        #Stores a specific grid point which is our x0 of the form (X[i][j], Y[i][j])
        fig, ax = plt.subplots()
        ax.pcolor(X,Y,A)

    def plot_iteration(self, N, a, b, c, d, simplified):
        x = np.linspace(a,b,N)
        y = np.linspace(c,d,N)
        X, Y = np.meshgrid(x,y)
        A = np.array([np.array([self.zeroes_index(np.array([X[i][j],Y[i][j]]), simplified)[1] for i in range (N)]) for j in range(N)])
        fig, ax = plt.subplots()
        ax.pcolor(X,Y,A)



class OG_class(fractal2D_erik, Fractal2D_amer):
     def __init__(self,F, J, f,j=None):
          fractal2D_erik.__init__(self, f, j=None)
          Fractal2D_amer.__init__(self, F, J)


     def choose_method(self, a,b,c,d,N, amer=None, erik=None):
          if amer and not erik:
               return self.Comlete(a,b,c,d,N)
          elif erik and not amer:
               return self.plot(N, a,b,c,d)








































