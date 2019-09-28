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
from sympy import symbols, diff

class Complex:
    def __init__(self, re, im, *args):
        if not re:
            self.re = lambda *args: self.R(*args)
        else:
            self.re = re


        if not im:
            self.im = lambda *args: self.I(*args)
        else:
            self.im = im


    def __repr__(self):
        if self.im < 0:
            return f"{self.re} - {self.im}i"
        elif self.im > 0:
            return f"{self.re} + {self.im}i"

    """-------------------------------------------------------------------------------------
    """

    def __add__(self, *other):
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        val1, val2 = self.re, self.im
        All = conversion_block(self,*other)
        for s in All:
            val1 += s.re
            val2 += s.im
        return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """

    def __sub__(self, *other):
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        val1, val2 = self.re, self.im
        All = conversion_block(self,*other)
        for s in All:
            val1 += s.re
            val2 += s.im
        return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """


    def __mult__(self, *other):
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        All_Complex = conversion_block(self, *other)
        start_self = [self]
        for s in All_Complex:
            new = Multiply(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]

    """------------------------------------------------------------------------------------
    """
    def __div__(self, *other, margin = 1.e-10):
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        def Divide(a,b):
            """
            this function divides thw instances of the class Complex
            """
            c = Multiply(a, -b)
            absb = float(b.re)**2 + float(b.im)**2
            return Complex(c.re/absb, c.im/absb)

        All_Complex = conversion_block(self, *other)
        start_self = [self]
        for s in All_Complex:
            abss = float(s.re)**2 + float(s.im)**2
            if abss < margin:
                raise Exception("div, other = {other} contains a number that is too small")
            new = Divide(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]

    """---------------------------------------------------------------------------------------
    """

    def __get_attr__(self, *other):
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def prel_attr(*a, Re=[], Im=[]):
            Re = [s.re for s in a]
            Im = [s.im for s in a]
            return [Re, Im]
        All_Complex = conversion_block(self, *other)
        return prel_attr(*All_Complex)

    """------------------------------------------------------------------------------------------
    """
    def __is_atrr__(self, *other):
        """
        THis metjod should give False or True depending on the input.
        (i am vey unsure if this method works properly)
        """
        def conversion_block(self, *other):
            """
            this function converts all the elements in other to instances
            of the class Complex
            """
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def ineq_(a, c=None):
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
                    return is_real_prel(a)
            elif np.abs(float(a.re)) > np.abs(float(a.im)):
                c == 2
                if c ==2:
                    return is_imag_prel(a)

        All_Complex = conversion_block(self, *other)
        empty = []
        for s in All_Complex:
            empty.append(ineq_(s))
        return empty




I = Complex(1,2)
P = Complex(2,3)
IP =

print(I)
































"""

class CN:
    def __init__(self, Re, Im):
        if not isinstance(Re, int) and not isinstance(Re, float):
            raise TypeError("the real part should be a number")
        if not isinstance(Im, int) and not isinstance(Im, float):
            raise TypeError("the imaginary part should be a number")
        self.Re = Re
        self.Im = Im

    def __repr__ (self):
        if self.Im < 0:
            return "{} - {}i".format(self.Re, -self.Im)
        return "{} + {}i".format(self.Re, self.Im)

    def __CNconj__(self):
        return CN(self.Re, -self.Im)

    def __RET_Re__(self):
        return np.array([self.Re, self.Re])[0]
    def __RET_Im__(self):
        return np.array([self.Im, self.Im])[-1]


    def __CNadd__(self, other):
        if type(other) == type(self.Re):
            return CN((self.Re+other) , self.Im)

        if type(other) == type(self.Im):
            return CN(self.Re , (self.Im+other))

        if type(other) == CN:
            return CN(self.Re + other.Re , self.Im + other.Im)


    def __CNmult__(self, other):
        a1, a2 = self.Re, self.Im
        b1, b2 = other.Re, other.Im
        return CN(a1*b1 - a2*b2, a1*b2 + a2*b1)

    def __CNdiv__(self, other):
        a1, a2 = self.Re, self.Im
        b1, b2 = other.Re, other.Im
        return CN((a1*b1 - a2*b2)/(b1**2 + b2**2), (a1*b2 + a2*b1)/(b1**2 + b2**2))


I = Complex(1,2)
print(I)

"""












"""


class Vectorspace:
	def __init__(self, basis):
		self.basis = np.array(basis)
#------------------------------------------------------
	def __add__(self, other):
		if type(other) == int:
            #det visar sig att denna funktion möjligtgör för en
            #att lägga till en integer till varje  (elementwise) i basen förvektorrummet.
            #SÅÅ detta betyder att dessa methods som hector gjorde
            #skapade verkar på själva bas(vetorerna, matriserna,...etc)
            #vi ser att additionen av en integer sker INUTI (self.basis + other)
			return type(self)(self.basis + other)

		if type(other) == Vectorspace:
			return type(self) (self.basis + other.basis)
#------------------------------------------------------
	def __radd__(self, other):
		if type(other) == int:
			return self + other



V = Vectorspace([[1,2,4],[0,4,1],[2,3,9]])
W = Vectorspace([[1,2,3], [5,6,8], [3,4,5]])
print((W+V + 5).basis)
print(V.__radd__(3))
print(Vectorspace.__radd__(V,3))




class Vektorrum:
    def __init__(self, bas):
        if not isinstance(bas, np.array()):
            raise TypeError("bas måste vara en vector")
        self.bas = np.array(bas)

#----ett vektorrum måste uppfylla kraven givna av axiomen för vektorrum
        #linjaritet
        #additivetet
        #skalning
        #identitet- och nollelementet
        #multiplicationsregel




class RationalNumber:
    def __init__(self, numerator, denominator):
        if not isinstance(numerator, int):
            raise TypeError("numerator should be of type in")
        if not isinstance(denominator, int):
            raise TypeError("denominator should be of type int")
        self.numerator = numerator
        self. denominator = denominator

    def convert2float(self):
        return float(self.numerator)/float(self.denominator)

    def add(self, other):
        p1, q1 = self.numerator, self.denominator
        if isinstance(other, RationalNumber):
            p2, q2 = other.numerator, other.denominator
        elif isinstance(other, int):
            p2, q2 = other, 1
        else:
            raise TypeError("...")
        return RationalNumber(p1*q2 + p2*q1, q1*q2)

        def __repr__ ( self ):
            return "{} / {}".format(self.numerator , self.denominator )

#------------------------
#Now we've created a method that initialises an instance with two attributes:
q = RationalNumber(10, 20)
#print("numerator =", q.numerator, "___", "denominator =", q.denominator)
#creating an adding method
#------------------------
#let's try the division method:
p = RationalNumber(10, 20)
S = p.convert2float()
print(S)
#------------------------
RationalNumber.convert2float(q)
#och
q.convert2float() # är ekvivalenta beteckningar
#------------------------
#let us try the adder method
q = RationalNumber(1, 2)
p = RationalNumber(1, 3)
c = RationalNumber.add(q,p)
print(c)




#import pandas as pd

#mu=1
#n=50
#dt=0.1
#x0=100
#x=pd.DataFrame()
#np.random.seed(1)

#for sigma in np.arange(0.8,2,0.2):
#    step=np.exp((mu-sigma**2/2)*dt)*np.exp(sigma*np.random.normal(0,np.sqrt(dt),(1,n)))
#    temp=pd.DataFrame(x0*step.cumprod())
#    x=pd.concat([x,temp],axis=1)

#x.columns=np.arange(0.8,2,0.2)
#plt.plot(x)
#plt.legend(x.columns)
#plt.xlabel('t')
#plt.ylabel('x')
#plt.title('Realizations of Geometric Brownian Motion with different variances\n mu=1')
#plt.show()

c









#-----------------bilda en method som multiplicerar två komplexa tal och
            #-----en methjod som dividerar två complexa tal

#----------testing zone-------------------
Z = CN(2,1)
E = CN.__RET_Re__(Z)
Q = CN(2,2)
D = CN.__CNadd__(Q,Q)
P = CN.__CNmult__(Z,Q)
Pt = CN.__CNconj__(P)
print(D)
print(E)
print(P)
print(Pt)
print(CN.__CNdiv__(Z,Q))

print(16*24)




def fac(x):
    if x == 0:
        return 1
    else:
        return x*fac(x-1)

def R(x,n):
    return (x-pi/3)**(2*n +1)/fac(2*n +1)

for k in range(5):
    if R(pi/3,k) > 1/2*10**(-3):
        print("n", "=", k , ":", R(pi/3,k))
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






class CTM:
    def __init__(self, X, B, C):
         if not isinstance(B, int) and not isinstance(B, float):
             raise TypeError("the real part should be a number")
         if not isinstance(C, int) and not isinstance(C, float):
            raise TypeError("the imaginary part should be a number")
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


"""













"""
def __add__(self, *other):
        def conversion_block(self, *other):
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        val1, val2 = self.re, self.im
        All = conversion_block(self,*other)
        for s in All:
            val1 += s.re
            val2 += s.im
        return Complex(val1, val2)

    def __sub__(self, *other):
        def conversion_block(self, *other):
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        val1, val2 = self.re, self.im
        All = conversion_block(self,*other)
        for s in All:
            val1 += s.re
            val2 += s.im
        return Complex(val1, val2)




    def __mult__(self, *other):
        def conversion_block(self, *other):
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def Multiply(a,b):
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        All_Complex = conversion_block(self, *other)
        start_self = [self]
        for s in All_Complex:
            new = Multiply(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]

    def __div__(self, *other, margin = 1.e-10):
        def conversion_block(self, *other):
            Number=[], py_C=[], self_C=[]
            Number = [Complex(s,0) for s in other if type(s)==float or type(s)==int]
            py_C   = [Complex(s.real, s.imag) for s in other if type(s)==complex]
            self_C = [s for s in other if type(s)==type(self)]
            All    = Number + py_C + self_C
            return All

        def Multiply(a,b):
                a1, a2 = a.re, a.im
                b1, b2 = b.re, a.im
                return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        def Divide(a,b):
            c = Multiply(a, -b)
            absb = float(b.re)**2 + float(b.im)**2
            return Complex(c.re/absb, c.im/absb)

        All_Complex = conversion_block(self, *other)
        start_self = [self]
        for s in All_Complex:
            abss = float(s.re)**2 + float(s.im)**2
            if abss < margin:
                raise Exception("div, other = {other} contains a number that is too small")
            new = Divide(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]
"""






















"""



    def __operation__(self, *other):
         def abs_(a):
             a1, a2 = float(a.re), float(a.im)
             A = np.sqrt(a1**2 + a2**2)
             return A
         def multiply(a,b):
             a1, a2 = a.re, a.im
             b1, b2 = b.re, b.im
             return Complex(a1*b1 - a2*b2, a1*b2 + a2*b1)
         def conversion_block(*other):
             for s in other:
                 if type(s)==complex:
                     s.re = s.real
                     s.im = s.imag
                 else:
                     raise Exception(f"other contains non-(float, int or complex) object")
             return other
         def __add__(*other):
             conv_other = conversion_block(*other)
             q1, q2 = self.re, self.im
             for s in conv_other:
                 if type(s)==float or type(s)==int:
                     s.re = s
                     s.im = 0
                 q1 += s.re
                 q2 += s.im
             return Complex(q1, q2)
         def __sub__(*other):
             conv_other = conversion_block(*other)
             q1, q2 = self.re, self.im
             for s in conv_other:
                 if type(s)==float or type(s)==int:
                     s.re = s
                     s.im = 0
                 q1 -= s.re
                 q2 -= s.im
             return Complex(q1, q2)
         def __mult__(*other):
             L =[]
             conv_other = conversion_block(*other)
             for s in conv_other:
                 if type(s)==float or type(s)==int:
                     l_new = float(self.re)*s
                     L.append(l_new)
                 elif isinstance(s,Complex):
                     l_new = multiply(self, s)
                     L.append(l_new)
             return L[-1]
         def __conj__(*other):
             if len(other) == 1 and isinstance(other[0], Complex):
                 a = other[0]
                 a1, a2 = float(a.re), float(a.im)
                 return Complex(a1, -a2)
             elif len(other) == 1 and type(other[0])==float or type(other[0])==int:
                 return other[0]

             elif len(other) > 1:
                 B = [Complex(s.re, -s.im) for s in other if isinstance(s, Complex)]
                 C = [s for s in other if type(s)==float or type(s)==int or type(s)==complex and s.imag ==0]
                 D = [Complex(s.real, -s.imag) for s in other if type(s)==complex]
                 return B+C+D


         def __abs__(*other):
             conv_other = conversion_block(*other)


             L = [abs_(s) for s in conv_other]
             return L





         def __div__(*other):
             A = __conj__(*other)
             L_div = [__mult__(self, __conj__(s)) for s in A]



             return



         return

"""


































































































































































