# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 18:48:50 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import symbols, diff
import time
from scipy.integrate import quad
from math import *

e = np.array([0., 1.,  2.,  33.])
print(e.shape)
E = np.array([[0., 1.,  2.,  33.],
              [1., 5.,  6.,  66.],
              [2., 9.,  10., 84.],
              [3., 12., 13., 99.],
              [4., 15., 16., 46.],
              [5., 18., 19., 87.],
              [6., 21., 22., 45.],
              [7., 24., 25., 77.],
              [8., 27., 28., 44.]])

e = np.array([1,1,1,1])

f = np.row_stack([e,E,e])
print(f)
f = list(f)
print(f)
print(type(f[3]))
print(type(f))
f = np.array(f)
print(type(f))





class Interval:
    def __init__(self, a,b):
        self.a = a
        self.b = b

        if not a:
            self.a = b
            self.b = b
        elif not b:
            self.a = a
            self.b = a

    def __repr__(self):
        return f"[{self.a}, {self.b}]"

    def __add__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.a + other[0].a, self.b + other[0].b)
        elif len(other) > 1:
            val1, val2 = self.a, self.b
            for s in other:
                val1 += s.a
                val2 += s.b
        return Complex(val1, val2)

    def __radd__(self, other):
        return self + other


    def __sub__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.a + other[0].b, self.b + other[0].a)
        elif len(other) > 1:
            val1, val2 = self.a, self.b
            for s in other:
                val1 -= s.b
                val2 -= s.a
        return Complex(val1, val2)



    def __mult__(self, *other):
        def Multiply(A,B):
            """
            This sub-function takes two datapoints and constructs
            two vectors out of their inerval bounds. THen it takes their
            tensor product. (https://en.wikipedia.org/wiki/Tensor_product)
            once that is done, it takes the min and  the max of the resulting matrix

            ON PARAMETERS:
                A (Interval)
                B (Interval)

            ON RETURN:
                epsilon (Interval)

            """
            A1,A2 = A.a, A.b
            B1, B2 = B.a, B.b
            alpha = np.array([A1, A2])
            beta = np.array([B1, B2])
            gamma = lapha * beta.reshape(-1,1)
            epsilon = Interval(gamma.min(), gamma.max())
            return epsilon

        if len(other) == 0:
            raise Exception("other should contain at least one interval")
        elif len(other) == 1:
            return Multiply(self, other[0])
        elif len(other) > 1:
            satrt_self =[self]
            for s in other:
                new = Multiply(start_self[-1], s)
                start_self.append(new)
            return start_self[-1]

    def __truediv__(self, *other, margin=1.e-10):
        def Divide(A,B):
            """
            This sub-function takes two datapoints and constructs
            two vectors out of their inerval bounds. THen it takes their
            tensor product (but instead of multiplying the vector elements,
            it divides the elements of the first vetor  with the seconds').
            (https://en.wikipedia.org/wiki/Tensor_product)
            once that is done, it takes the min and  the max of the resulting matrix

            ON PARAMETERS:
                A (Interval)
                B (Interval)

            ON RETURN:
                epsilon (interval)

            """
            A1,A2 = A.a, A.b
            B1, B2 = B.a, B.b
            alpha = np.array([A1, A2])
            beta = np.array([B1, B2])

            if (beta==0).any():
                raise Exception("function Divide, can't divide by zero")

            else:
                gamma = lapha * beta.reshape(-1,1)
                epsilon = Interval(gamma.min(), gamma.max())
                return epsilon

        if len(other) == 0:
            raise Exception("other should contain at least one interval")
        elif len(other) == 1:
            return Divide(self, other[0])
        elif len(other) > 1:
            start_self =[self]
            for s in other:

                if np.abs(float(s.a)) < margin:
                    raise Exception(f"{s.a} is too small")
                elif np.abs(float(s.b)) < margin:
                    raise Exception(f"{s.b} is too small")
                elif float(s.a) < 0 and float(s.b) >0:
                    raise Exception(f"[{s.a}, {s.b}] contains zero")


                else:
                    new = Divide(start_self[-1], s)
                    start_self.append(new)
            return start_self[-1]



        def __pow__(self, n):

            L = [float(self.a)**n, float(self.b)**n]
            if not isinstance(n, int):
                raise Exception("n should be an ineger")

                def Pow(x,n, c="normal"):
                    x_an, x_bn = float(x.a)**n, float(x.b)**n
                    if c == "normal":
                        return Interval(x_an, x_bn)
                    elif c == "inv":
                        return Interval(x_bn, x_an)

            if n >= 1 and n%2 != 0:
                return Pow(self, n)

            elif n>0  and n%2 ==0:
                if a >= 0:
                    return Pow(self, n)
                elif b < 0:
                    return Pow(self, n, c="inv")
                else:
                    return Interval(0,max(L))



        def __Contains(self, other):
            def Con(S, M):
                a,b = S.a, S.b
                p,q = M.a, M.b

                if p>= a and q<=b:
                    return True
                elif p<a and q>b:
                    return f"inversely, S={S} contains M={M}"
                elif p<a and q==a:
                    return f"S={S} and M={M} intersect at 1 point; {S.a}"
                elif p==b and q>b:
                    return f"S={S} and M={M} intersect at 1 point; {S.b}"
                elif p<b and q>b:
                    return f"partly, their intersection is the interval {Interval(p,b)}"
                elif p<a and q>a:
                    return f"partly, their intersection is the interval {Interval(a,q)}"
                else:
                    return False


            if isinstance(other, int) or isinstance(other, float):
                l,u = self.a, self.b
                if l<= other <= u:
                    return True
                else:
                    return False


            if isinstance(other, list) and len(other)==2 :
                other = Interval(other[0], other[1])
                return Con(self, other)
            elif len(other) != 2:
                raise Exception("{other} cannot be turned into Interval type")



            if isinstance(other, np.ndarray) and other.shape[0]==2:
                other = Interval(other[0], other[1])
                return Con(self, other)
            elif other.shape[0] != 2:
                raise Exception("{other} cannot be turned into Interval type")



            if isinstance(other, tuple) and len(other)==2:
                other = Interval(other[0], other[1])
                return Con(self, other)
            elif len(other) != 2:
                raise Exception("{other} cannot be turned into Interval type")

            if isinstance(other, Interval):
                return Con(self,other)



        def __Plot__(l,u, N, delta, p=None):
            if not p:
                p = lambda I: 3*I**3 -2*I**2 - 5*I -1

            def PLOT(l, u, N, delta):
                vp = np.vectorize(p)
                xl = np.linspace(l,u,N)
                xu = np.linspace(l,u,N) + delta
                Tl = vp(xl)
                Tu = vp(xu)

                plt.figure()
                plt.plot(xl, Tl, label = "transformed lower bounds")
                plt.plot(xl, Tu, label = "transformed upper bouds")
                plt.xlabel("x")
                plt.ylabel("p(I)")
                plt.title("$p(I) = {p}$, I=interval(x,x+{delta})")
                plt.xlim(0,1)
                plt.ylim(-10,4)

                plt.title(f"p(I) = {p} = interval(x, x+{delyta})")
                plt.legend(loc ="best")

            if not l and not u and not N and not delta:
                l = 0
                u = 1
                N = 1000
                delta = 0.5
                PLOT(l,u,N,delta)
            else:
                PLOT(l,u,N,delta)














































class Complex:
    def __init__(self, R, I,x,y, re=None, im=None):
            self.re = R(x,y)
            self.im = I(x,y)

    def __repr__(self,re=None, im=None):
        return f"{self.re} + {self.im}i"


    def __add__(self, other, x,y, re=None, im=None, ):
        def __addd__(A,B):
            A.re = lambda x,y: A1(x,y)
            A.im = lambda x,y: A2(x,y)
            B.re = lambda u,v: B1(u,v)
            B.im = lambda u,v: B2(u,v)
            C1 = lambda x,y: A1(x,y) + B1(x,y)
            C2 = lambda x,y: A2(x,y) + B2(x,y)


            return Complex(C1, C2, x,y)

        return add(self, other)

A = Complex(lambda x,y: x, lambda x,y: y,1,2)
B = Complex(lambda u,v: u, lambda u,v: v,2,3)

print(A)
print(Complex.__add__(A,B,1,2,2,3))











    def __repr__(self):
        if self.im < 0:
            return f"{self.re} - {(-1)*self.im}i"
        else:
            return f"{self.re} + {self.im}i"
    """-------------------------------------------------------------------------------------
    """


    def __add__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.re + other[0].re, self.im + other[0].im)
        elif len(other) > 1:
            val1, val2 = self.re, self.im
            for s in other:
                val1 += s.re
                val2 += s.im
        return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """
    def __sub__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.re - other[0].re, self.im - other[0].im)
        elif len(other) > 1:

            val1, val2 = self.re, self.im
            for s in other:
                val1 -= s.re
                val2 -= s.im
            return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """
    def __mul__(self, *other):
        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        if len(other) == 0 :
            return 0
        elif len(other) == 1:
            return Multiply(self, other[0])
        elif len(other) > 1:
            start_self = [self]
            for s in other:
                new = Multiply(start_self[-1], s)
                start_self.append(new)
                return start_self[-1]

    """------------------------------------------------------------------------------------
    """
    def __div__(self, *other, margin = 1.e-10):
        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)
        def conj(a):
            a1, a2 = a.re, a.im
            return Complex(a1, -a2)

        def Divide(a,b):
            """
            this function divides thw instances of the class Complex
            """

            b_conj = conj(b)
            c = Multiply(a, b_conj)
            absb = float(b_conj.re)**2 + float(b_conj.im)**2
            return Complex(c.re/absb, c.im/absb)

        start_self = [self]
        for s in other:
            abss = float(s.re)**2 + float(s.im)**2
            if abss < margin:
                raise Exception("div, other = {other} contains a number that is too small")
            new = Divide(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]

    """---------------------------------------------------------------------------------------
    """

    def __get_attr__(self, *other):
        def prel_attr(self,*other, Re=[], Im=[]):
            Re = [s.re for s in other]
            Im = [s.im for s in other]
            return [Re, Im]

        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return prel_attr(self, other[0])
        elif len(other) > 1:
            return prel_attr(self,*other)

    """------------------------------------------------------------------------------------------
    """
    def __is_attr__(self, *other):
        """
        I do't think this method is very good. Needs revision.
        """
        def ineq_(a, c):
            def is_real_prel(a):
                if a.im == 0:
                    return True
                else:
                    return None
            def is_imag_prel(a):
                if a.re == 0:
                    return True
                else:
                    return None
            if np.abs(float(a.re)) < np.abs(float(a.im)):
                c ==1
                if c==1:
                    return is_real_prel(a,c)
            elif np.abs(float(a.re)) > np.abs(float(a.im)):
                c == 2
                if c ==2:
                    return is_imag_prel(a,c)

        empty = []
        for s in other:
            empty.append(ineq_(s,c))
        return empty

"""
f = lambda x,y: x
g = lambda x,y: y

A = Complex(lambda x,y: x, lambda x,y: y,  1,1)
print(A)
B = Complex(lambda x,y: x, lambda x,y: y,  2,2)
List = [B]
a = Complex.__add__(A, B)
"""



"""

S = Complex(3,3)
Q = Complex(1,2)
P = Complex(2,3)
other = [Q, P, S]
a = Complex.__add__(*other)
b =Complex.__sub__(*other)
c =Complex.__mul__(*other)
d = Complex.__div__(*other)
e = Complex.__get_attr__(*other)
f = Complex.__is_attr__(*other)
print([a,b,c,d,e,f])

"""
































