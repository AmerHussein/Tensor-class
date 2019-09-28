# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 09:23:16 2019

@author: Amer
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
import matplotlib.widgets as mpw
from matplotlib.widgets import Slider
from scipy.integrate import quad
import scipy.linalg as sl
from sympy import sqrt, sin, cos

help(enumerate)

def level_sets(f,level_values, x1,x2, y1,y2, N, tol=1.e-2):
     X=list(np.linspace(x1,x2, N))
     Y=list(np.linspace(y1,y2, N))











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


A = np.zeros((9,4))
E = np.array([[0., 1.,  2.,  33.],
              [1., 5.,  6.,  66.],
              [2., 9.,  10., 84.],
              [3., 12., 13., 99.],
              [4., 15., 16., 46.],
              [5., 18., 19., 87.],
              [6., 21., 22., 45.],
              [7., 24., 25., 77.],
              [8., 27., 28., 44.]])
print(E.reshape(4,9))
print(E[:,1])
A = A.reshape(4,9)
A.reshape(4,9)[0] = E[:,1]
print(A)

A = np.array([1,2,3])
B = np.zeros(0)
print(np.hstack([A, B]))
print(A.shape)





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
            return np.linalg.det(M),"linearly dependent 1"

        if np.linalg.norm(conclusion) < margin:
            return f"det = {np.linalg.det(M)}",conclusion, "linearly independent 2"
        else:
            return M,conclusion, "linearly dependent 3"

S = [np.array([1,1,1.2]), np.array([3,2.2,5])]
print(lin_dep(S))













a = [1]
print(a[0:1])
print(a[:0] + a[-0:])




def missing_numbers(a):
    E =[]
    a = list(a)
    while len(a) != 2:
        ind = a.index(max(a))
        E.append(max(a))
        a = a[:ind] + a[-(len(a) -(ind+1)):]

    E.append(max(a))
    E.append(min(a))
    E = E[::-1]
    E_min = min(E)

    if E_min < 0:
        raise Exception("no negative numbers allowed")
    else:
        b = np.arange(E_min)
        b = list(b)
    i=0
    while i < len(E) -1:
        if E[i+1] -E[i] > 1:
            for j in range(1, E[i+1] -E[i]-1):
                b.append(E[i] +j)
            i+=1
        elif E[i+1] -E[i] ==1:
            i+=1
        elif E[i+1] -E[i] ==0:
            E = E[:i] + E[-(len(E) -(i+1)):]
            i+=1
    return b

print(missing_numbers([4, 5, 8, 9]))




def num_coinse(a, B, object=None):
    """
    This function gives you back the right amount of change.
    You can choose your own denominations regarding the type of coins
    available in a chash register.

    ON INPUTS:
        a (integer):
            the change
        B (list):
            the denominations
    ON RETURN:
        (e, L) (int, list):
            the amount of coins and a list containing information
            about which type coins were used and how many of each
    """
    if B[0] == object:
        B = B[-(len(B)-1):]
    if B[-1] == object:
        raise Exception("nigga, you need to have 1")
    if object in B:
        E = [B,[]]
        i=0
        while object in E[0]:
            if E[0][0] == object:
                E[0] = E[0][-(len(E[0])-1):]
            elif E[0][i]==object:
                E[1].append(E[0][:i])
                E[0] = E[0][-(len(E[0])- (i+1)):]
            else:
                i+=1
        E[1].append(E[0])


        for i in range(1, len(E[1])):
            E[1][0] += E[1][i]
        B = E[1][0]

    L = [B,[[0] for i in range(len(B))]]
    k=0
    while a > 0:
        if a>=B[k]:
            a-=B[k]
            L[1][k][0]+=1
        else:
            k+=1
    e = 0
    for i in range(len(B)):
           e += L[1][i][0]
    return e, L

print(num_coinse(76, [None,50, 30,None, 20, 10, None,None, None, None,None, None, None, None,None, None, None, None,None, None, None, None,None, None, None, None,None, None, None, None,None, None,7, 5, None, 2, 1]))



def PLOT(a,b):
    " run the code and search ""mandelbrot set"" on wikipedia"
    def bifurcation(r,x0 ,max_iter=1000):
        L =[x0]
        k=0
        while k < max_iter:
            new = r*L[-1]*(1-L[-1])
            L.append(new)
            k+=1
        return L[-1]

    fig, ax = plt.subplots()
    R = np.linspace(1,4, a)
    X = np.linspace(0,1, b)
    for x0 in X:
        ax.plot(R, bifurcation(R, x0), label = "bifurvation for r")

print(PLOT(100, 100))










fig = figure ()
sld_ax = axes([0.2, 0.9, 0.6, 0.05])
plt_ax = axes([0.1, 0.15 , 0.8, 0.7])
val_init=1.0
sld = Slider(sld_ax , "values", 0., 5.,valinit=val_init ,   valfmt="%1.1f")
anno=plt_ax.annotate(f"The  given  value is{val_init:1.1f}",(0.2,0.40),fontsize=24)
def val_annotate(val):
    anno.set_text(f"The  given  value is {val:1.1f}")
sld.on_changed(val_annotate)

fig = plt.figure(figsize = (4,2))
sld_ax = plt.axes([0.2, 0.9, 0.6, 0.05])# axes  for  slider
ax = plt.axes([0.1, 0.15 , 0.8, 0.7])# axes  for the  plot
sld = Slider(sld_ax , "amp", 0., 5.)
x = np.linspace(-2*np.pi, 2*np.pi, 200)
ax.set_ylim(-5.5, 5.5)
ax.plot(x, np.cos(x))
lines, = ax.plot(x, 1-((x-sld.val)**2)/2)
def update_amplitude(val):
    lines.set_ydata(val)

sld.on_changed(update_amplitude)



"""
GUI: Graphical User Interface: using matplotlib
    you could change the data attribute to change the influnce of parameters of our curves.
    "maybe making a slider animation"
    - a widget is a little window which serves as aa user interface to a computer program.

    we make a slider.
    1) make a figure object fig = figure()
    2) make an axes object wthin the figure ax =axes(<list type object>)
    3) Create a widget, e.g. a slider, by using this axes object sld = Slider(ax,...)
    the signature of the axes-command:
    ax = axes([left, nottom, width, heuight]) - this is more of an aesthetic
    4) we turn this rectangle that is created by axes([<>])
    sld = Slider(ax, "attribute", float(max), float(min), valinit=1, valfmt="%1.1"f)





"""





"""
class FlatFunctionError(Exception):
    pass
    """"""
    the class FlatFunctionError , you don't need to add new methods
    to the derived class to be able to rasie completelly new Errors.
    WHat you do is tat when you are aware that you code will encounter
    an error that isn't built-in in pythin, you simply write, as an example
    """"""
class NUMA01Exception(Exception):
    pass

def f1(x):
    c = f2(x)
    return c
def f2(x):
    try:
        c = f3(x)
    except ZeroDivisionError:
        raise NUMA01Exception("this is a testfunction")
    return c
def f3(x):
    return 1/x

f1(0)
"""

class Interval:
    def __init__(self, *args):
        """
        args should be a tuple. and if you want a tuple with only one
        element then you could write args = (a, ) where a is a number.
        """
        if len(args) == 1:
            self.a = args[0]
            self.b = args[0]

        elif len(args) == 2:
            self.a = args[0]
            self.b = args[1]

    def __repr__(self):
        return f"[{self.a}, {self.b}]"

    def __add__(self, *other):
        v1, v2 = self.a, self.b

        if len(other)==0:
            return self

        elif len(other)==1:
            G = other[0]
            if isinstance(G, int) or isinstance(G, float):
                s = Interval(0, G)
                return Interval(v1 + s.a, v2 + s.b)

            elif isinstance(G, Interval):
                return Interval(v1 + G.a, v2 + G.b)

        elif len(other) > 1:
            for s in other:
                if isinstance(s, float) or isinstance(s, int):
                    s = Interval(0, s)
                v1 += s.a
                v2 += s.b

            return Interval(v1, v2)

    def __radd__(self, *other):
        v1, v2 = self.a, self.b

        if len(other)==0:
            return self

        elif len(other)==1:
            G = other[0]
            if isinstance(G, int) or isinstance(G, float):
                s = Interval(G, 0)
                return Interval(v1 + s.a, v2 + s.b)

            elif isinstance(G, Interval):
                return Interval(v1 + G.a, v2 + G.b)

        elif len(other) > 1:
            for s in other:
                if isinstance(s, float) or isinstance(s, int):
                    s = Interval(s,0)
                v1 += s.a
                v2 += s.b

            return Interval(v1, v2)




    def __sub__(self, *other):
        v1, v2 = self.a, self.b

        if len(other)==0:
            return self

        elif len(other)==1:
            G = other[0]
            if isinstance(G, int) or isinstance(G, float):
                s = Interval(0, G)
                return Interval(v1 - s.b, v2 - s.a)

            elif isinstance(G, Interval):
                return Interval(v1 - G.b, v2 - G.a)

        elif len(other) > 1:
            for s in other:
                if isinstance(s, float) or isinstance(s, int):
                    s = Interval(0, s)
                v1 -= s.b
                v2 -= s.a

            return Interval(v1, v2)

    def __rsub__(self, *other):
        v1, v2 = self.a, self.b

        if len(other)==0:
            return self

        elif len(other)==1:
            G = other[0]
            if isinstance(G, int) or isinstance(G, float):
                s = Interval(G,0)
                return Interval(v1 - s.b, v2 - s.a)

            elif isinstance(G, Interval):
                return Interval(v1 - G.b, v2 - G.a)

        elif len(other) > 1:
            for s in other:
                if isinstance(s, float) or isinstance(s, int):
                    s = Interval(s,0)
                v1 -= s.b
                v2 -= s.a

            return Interval(v1, v2)


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
            gamma = alpha * beta.reshape(-1,1)
            epsilon = Interval(gamma.min(), gamma.max())
            return epsilon

        if len(other) == 0:
            raise Exception("other should contain at least one interval")
        elif len(other) == 1:
            return Multiply(self, other[0])
        elif len(other) > 1:
            start_self =[self]
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
                raise ZeroDivisionError("child-function; Divide: can't divide by zero")
            elif beta[0] < 0 and beta[1] > 0:
                raise Exception("the _nämnaren_ contains zero")
            else:
                gamma = alpha / beta.reshape(-1,1)
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
            raise Exception("n should be an integer")

        def Pow(x,n, c="normal"):
            x_an, x_bn = float(x.a)**n, float(x.b)**n
            if c == "normal":
                return Interval(x_an, x_bn)
            elif c == "inv":
                return Interval(x_bn, x_an)

        if n >= 1 and n%2 != 0:
            return Pow(self, n)

        elif n>0  and n%2 ==0:
            if self.a >= 0:
                return Pow(self, n)
            elif self.b < 0:
                return Pow(self, n, c="inv")
            else:
                return Interval(0,max(L))



    def __contain__(self, other):
        def Con(S, M):
            a,b = S.a, S.b
            p,q = M.a, M.b
            if p>= a and q<=b:
                return True
            elif p<a and q>b:
                return f"inversely, {S} contains {M}"

            elif p<a and q==a:
                return f"{S} and {M} intersect at the point x={S.a}"
            elif p==b and q>b:
                return f"{S} and {M} intersect at the point x={S.b}"
            elif a<p<b and q>b:
                return f"partly, their intersection is the interval {Interval(p,b)}"
            elif p<a and b>q>a:
                return f"partly, their intersection is the interval {Interval(a,q)}"
            else:
                return False


        if isinstance(other, int) or isinstance(other, float):
            l,u = self.a, self.b
            if l<= other <= u:
                return True
            else:
                return False


        elif isinstance(other, np.ndarray):
            if other.shape[0] != 2:
                raise Exception("{other} cannot be turned into Interval type")
            else:
                other = Interval(other[0], other[1])
                return Con(self, other)


        elif isinstance(other, tuple) or isinstance(other, list):
            if len(other) !=2:
                raise Exception("{other} cannot be turned into Interval type")
            else:
                other = Interval(other[0], other[1])
                return Con(self, other)


        elif isinstance(other, Interval):
            return Con(self,other)



def Plot(l,u, N, delta, p=None):
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
        plt.plot(xl, Tu, label = "transformed upper bounds")
        plt.xlabel("x")
        plt.ylabel("p(I)")
        plt.title("$p(I) = {p}$, I=interval(x,x+{delta})")
        plt.xlim(0,1)
        plt.ylim(-10,4)

        plt.title(f"p(I) = {p} = interval(x, x+{delta})")
        plt.legend(loc ="best")

    if not l and not u and not N and not delta:
        l = 0
        u = 1
        N = 1000
        delta = 0.5
        PLOT(l,u,N,delta)
    else:
        PLOT(l,u,N,delta)


E = Interval(1,)
print(E)
A = Interval(8,5)
Q = Interval(5,2)
q = Interval(1,1)
print(E-4)
print(Q)
print(A)
L =[q,q,q,q,q,q,q,q,q,q]
print(Interval.__add__(Q,q))
print(Interval.__truediv__(Q,A))
print(A-q)
Plot(l=None, u=None, N=None, delta=None)





























































































































































"""






class interval:
    def __init__(self, *args):
        r = args[0]
        if len(args)>1:
            z = args[1]
            self.l_b = r
            self.u_b = z
        elif len(args)==1:
            self.l_b = r
            self.u_b = r

    def __repr__ (self):
        M = [self.l_b, self.u_b]
        return "[{} , {}]".format(M[0], M[1])

    def __add__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]+N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(M[0]+sum(N))
                Q.append(M[1]+sum(N))
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)

    def __radd__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]+N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(sum(N)+M[0])
                Q.append(sum(N)+M[1])
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)
    def __sub__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]-N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(M[0]-sum(N))
                Q.append(M[1]-sum(N))
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)

    def __rsub__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]-N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(sum(N)-M[0])
                Q.append(sum(N)-M[1])
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)
    def __mul__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]*N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(M[0]*sum(N))
                Q.append(M[1]*sum(N))
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)


    def __rmul__(self, other):
        M = np.array([self.l_b, self.u_b])
        Q = []
        if type(other)==Interval:
            self.other = np.array([other.l_b, other.u_b])
            self.other = list(self.other)
            N = self.other
            for i in range(2):
                for j in range(2):
                    Q.append(M[i]*N[j])
            l_b = min(Q)
            u_b = max(Q)
        if type(other)==int or float:
            if type(other)!=Interval:
                r = [other]
                N = [r]
                Q.append(sum(N)*M[0])
                Q.append(sum(N)*M[1])
                l_b = min(Q)
                u_b = max(Q)
        return Interval(l_b, u_b)
    def __truediv__(self, other, tol=10**(-4)):
        M = np.array([self.l_b, self.u_b])
        N = np.array([other.l_b, other.u_b])
        R = []
        for i in range(2):
            if other.l_b == 0:
                print("Division with zero!")
                break
            elif other.l_b <1.e-50 and other.l_b >-1.e-50:
                print("Resulting interval is infinitely large!")
                break
            elif other.u_b ==0:
                print("Division with zero!")
                break
            elif other.u_b <1.e-50 and other.u_b >-1.e-50:
                print("Resulting interval is infinitely large!")
                break
            else:
                for j in range(2):
                    R.append(M[i]/N[j])
                    l_b = min(R)
                    u_b = max(R)
        return Interval(l_b, u_b)

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

    def __pow__(self, n):
        M = np.array([self.l_b, self.u_b])
        if M[0] >= 0:
            c = [M[0]**n, M[1]**n]
            return Interval(*c)
        elif M[1] < 0:
            c = [M[1]**n, M[0]**n]
            return Interval(*c)
        else:
            c = [M[0]**n, M[1]**n]
            d = [0, max(c)]
            return Interval(*d)

def createintervals(self, xl, xu):
        I=[]
        for i in range(len(xl)):
            l_b=xl[i]
            u_b=xu[i]
            I.append(Interval(l_b, u_b))
        return I

def polynomial(self, I):
    p=[]
    for i in range(1000):
        polynomial=(((3*(I[i]**3))-(2*(I[i]**2))-(5*(I[i])))-1)
        p.append(polynomial)
    return p


xl=np.linspace(0.,1,1000)
xu=np.linspace(0.,1,1000) + 0.5
tempInterval=Interval(0, 0)
I=tempInterval.createintervals(xl, xu)
p=tempInterval.polynomial(I)

yl=[]
yu=[]
for i in range(1000):
    yl.append(p[i].l_b)
    yu.append(p[i].u_b)


plt.plot(xl, yl)
plt.plot(xl, yu)
plt.xlabel("x")
plt.ylabel("p(I)")
plt.title("$p(I) = 3I^3 − 2I^2 − 5I − 1$, I=interval(x,x+0.5)")
plt.xlim(0,1)
plt.ylim(-10,4)
plt.show()

"""

















































