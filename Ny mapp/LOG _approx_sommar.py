# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 11:45:28 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
from scipy.integrate import quad,dblquad
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import scipy.linalg as sl

def log_approx(x,n,difference = "off"):
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

    if difference == "off":
        return (x-1)/ a[-1]
    elif difference == "on":
        return np.abs((x-1)/a[-1] -np.log(x))



def fastlogapprox_amer(x,n, fast_diff_mat="off"):

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
            d[j,i] = (d[j-1,i] -(d[j-1,i-1]*4**(-j)))/((1-4**(-j)))

    if fast_diff_mat == "on":
        d[-1:,0] = np.abs((x-1)/d[-1,-1] - (x-1)/a[-1])
        return d
    elif fast_diff_mat=="off":
        return (x-1)/d[-1,-1]

VALUES_fast =[]
VALUES_slow =[]
DIFFERENCE = []
X = np.linspace(2,100,100)
for i in range(10):
    VAL1 = [fastlogapprox_amer(x,i) for x in X]
    VAL2 = [log_approx(x,i) for x in X]
    diff = [np.abs(fastlogapprox_amer(x,i) -log_approx(x,i)) for x in X]
    DIFFERENCE.append(diff)
    VALUES_fast.append(VAL1)
    VALUES_slow.append(VAL2)

plt.figure()
for s in DIFFERENCE:
    plt.plot(X, s, label = f"{i} iterations")
plt.legend(loc = "best")
plt.title("difference in the value of the fast/sow functions different n")


plt.figure()
for s in VALUES_fast:
    plt.plot(X, s, label = f"{i} iterations")
plt.legend(loc = "best")
plt.title("plotting the fast function for different n")

plt.figure()
for s in VALUES_slow:
    plt.plot(X, s, label = f"{i} iterations")
plt.legend(loc = "best")
plt.title("plotting the slow function for different n")

E1 = [fastlogapprox_amer(5,i) for i in range(10)]
E2 = [log_approx(5,i) - np.log(5) for i in range(10)]
E3 = []
e = [i for i in range(10)]

plt.plot(e, E2, label=  "error of the slow approximation")
plt.figure()
plt.plot(e, E1, label ="error of charlson acceleerated method")


plt.title("comparing the error of the fast and slow approximation")
plt.legend(loc = "best")







plt.figure()
for i in range(4):
    plt.plot(X, fastlogapprox_amer(X,i), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the function for different n")










def fastlogapprox_hector(x, n):
    a = np.zeros(n)
    g = np.zeros(n)
    d = np.zeros((n,n))
    a[0] = 0.5*(1+x)
    g[0] = sqrt(x)
    for i in range(1,n):
        a[i] = 0.5*(a[i-1]+g[i-1])
        g[i] = np.sqrt(a[i]*g[i-1])
    d[0,:] = a
    for i in range(0,n):
        d[i,i] = (a[i] - 4**(-i)*a[i-1])/(1-4**(-i))
    return d,(x-1)/d[-1,-1]



def fastlogapprox_alex(x, n):
    a0 = 0.5*(1+x)
    g0 = np.sqrt(x)
    a = [a0]
    g = [g0]
    d = np.zeros((n,n))
    for i in range(n):
        if i<n:
            a.append(0.5*(a[i]+g[i]))
            g.append(np.sqrt(a[i+1]*g[i]))
            d[0,i] = a[i]
        elif i==n:
            break
    for k in range(1, n):
        for i in range(1, n):
            d[k,i] = ((d[k-1,i] - 4**(-k)*d[k-1, i-1])/(1-4**(-k)))
        return d,(x-1)/d[k,i]





A = fastlogapprox_amer(5,5)

print(A[-1,-1])

print(fastlogapprox_hector(5, 5))
print("-------------------------------------------------------------")
print(fastlogapprox_alex(5, 5))
print("-------------------------------------------------------------")
print(fastlogapprox_amer(5, 4))



def power_decomposition(a):
        """
        This function takes in a number a and returns the part of the number
        that can be decomposed into a linear combination of powers of 10.

        ON PARAMETERS:
            a (float)

        ON RETURN:
            power_decomp (list):
                this list contains elements whose index values are the same as the
                exponent value that 10 is suposed to raised to.

                EXAMPLE:
                    a = 1334.976
                    resdiue(a)
                    >>> [4,3,3,1]

        ATTENTION!:
            This function does not take into account the decimals of the input. (this was intentional though)

            THIS FUNCTION IS VERY BAD WHEN DEALING WITH LARGE NUMBERS AND NUMBERS
            THAT END WITH MORE THAN TWO ZEROS FOR SOME REASON , IT HAS SOME MAJOR FLAWS AND BUGS THAT
            MAKE THE ANSWER RESULT WEIRD (THE AMOUNT OF DIGITS SUDDENLY BECOME MORE THAN THE AMOUNT OF
            DIGITS OF THE INPUT AND SOME DIGITS MAY BE REPLACED BY "RADNOM" DIGITS) .

            BUT IT WORKS FOR SMALLER NUMBERS WITH LESS THAN 16 DIGITS .

        """
#        K = np.array([10**(-j) for j in range(40)])
        L = np.array([10**(j) for j in range(40)])

        power_index  = np.where(a > L)
        power_value = len(list(power_index[0]))
        power_decomp =[]
        absa = np.abs(a)
        for i in range(1,power_value+1):
            for j in range(1000):
                absa -= 10**(power_value-i)
                if absa < 0:
                    power_decomp.append(j)
                    absa = absa + 10**(power_value-i)
                    break
                if absa == 0:
                    power_decomp.append(0)
                if  1<= np.abs(absa) <= 2:
                    power_decomp.append(1)
        if power_decomp[0] == 0:
            power_decomp[0] = 1
        return power_decomp[::-1]
print(power_decomposition(863))





def rounding_(a, decimal_place, c ="one"):
    """
    This function lets you round a number (a) to a given decimal_place

    ON INPUT:
        a (float)

    """
    def power_decomposition(a):
        """
        This function takes in a number a and returns the part of the number
        that can be decomposed into a linear combination of powers of 10.

        ON PARAMETERS:
            a (float)

        ON RETURN:
            power_decomp (list):
                this list contains elements whose index values are the same as the
                exponent value that 10 is suposed to raised to.

                EXAMPLE:
                    a = 1334.976
                    resdiue(a)
                    >>> [4,3,3,1]
        ATTENTION!:
            This function does not take into account the decimals of the input. (this was intentional though)

            THIS FUNCTION IS VERY BAD WHEN DEALING WITH LARGE NUMBERS AND NUMBERS
            THAT END WITH MORE THAN TWO ZEROS FOR SOME REASON , IT HAS SOME MAJOR FLAWS AND BUGS THAT
            MAKE THE ANSWER RESULT WEIRD (THE AMOUNT OF DIGITS SUDDENLY BECOME MORE THAN THE AMOUNT OF
            DIGITS OF THE INPUT AND SOME DIGITS MAY BE REPLACED BY "RADNOM" DIGITS) .

            BUT IT WORKS FOR SMALLER NUMBERS WITH LESS THAN 16 DIGITS .

        """
#        K = np.array([10**(-j) for j in range(40)])
        L = np.array([10**(j) for j in range(40)])

        power_index  = np.where(a > L)
        power_value = len(list(power_index[0]))
        power_decomp =[]
        absa = np.abs(a)
        for i in range(1,power_value+1):
            for j in range(1000):
                absa -= 10**(power_value-i)
                if absa < 0:
                    power_decomp.append(j)
                    absa = absa + 10**(power_value-i)
                    break
                if absa == 0:
                    power_decomp.append(0)
#                if  1<= np.abs(absa) <= 2:
 #                   power_decomp.append(1)
 #       if power_decomp[0] == 0:
  #          power_decomp[0] = 1
        return power_decomp[::-1]

    def reverse_decomposition(A):
        """
        this function takes a list and uses its elements (numbers)
        and uses them as wheights to be multiplied with 10**(i) where i is the
        element's index inside the list.
        (it basically reverses the procedure of the function called: "power_decomposition",
        and returns a but with the decimals subtrracted away)


        ON INPUT:
            A (list)

        ON RETURN:
            val (int)

        EXAMPLE:
             A = [4,3,3,1]
             reverse_decomposition(A)
             >>> 1334


        """
        if type(A) != list:
            raise Exception("The object to be reverse decomposed should be a list")
        l = len(A)
        V = A
        val = 0
        for i in range(l):
            val += V[i]*10**(i)
        return val

    pd                    = power_decomposition(a)
    a_without_decimals    = reverse_decomposition(pd)
    decimals              = a - a_without_decimals
    decimals_mult         = decimals * 10**(decimal_place+1)
    decimals_LIST         = power_decomposition(decimals_mult)
    dm                    = decimals_LIST

    if len(dm) ==1:
#        if dm[0] >= 5:
 #           dm[0] = 0
         if decimal_place == 0:
             if decimals_mult >= 5:
                 return reverse_decomposition(power_decomposition(a+1)), decimals_mult
             elif decimals_mult < 5:
                 return  reverse_decomposition(power_decomposition(a))

    elif len(dm) == 0:
        if decimal_place == 0:
                if decimals_mult >= 5:
                    return reverse_decomposition(power_decomposition(a+1)), decimals_mult
                elif decimals_mult < 5:
                    return  reverse_decomposition(power_decomposition(a))



    elif len(dm) != 1:
        if dm[0] >= 5:
            dm[0] = 0
            dm[1] = dm[1] +1

        elif dm[0] < 5:
            dm[0] = 0
        decimals_rounded_mult = reverse_decomposition(dm)
        decimals_rounded      = decimals_rounded_mult / 10**(decimal_place+1)
        rounded_a             = decimals_rounded + a_without_decimals

        if c == "many":
            return f"{rounded_a}{0} +- {10**(-decimal_place)/2}" , rounded_a, a, dm
        elif c =="one":
            return rounded_a, power_decomposition(a), len(power_decomposition(a))

print(rounding_(3.4376540, 3,c="one"))





def derivative_func(f=None, *args, h= 1.e-10):
    if not f:
        f = lambda *args: func(*args)

    def finite_diff(f, s):
        deriv_f = (f(s + h) - f(s))/h
        return deriv_f

    if type(f) != np.array():
        Grad_f = np.array([finite_diff(f, s) for s in args])
        return Grad_f









def power_decomposition(a):
    """
    This function takes in a number a and returns the part of the number
    that can be decomposed into a linear combination of powers of 10.

    ON PARAMETERS:
        a (float)

    ON RETURN:
        power_decomp (list):
            this list contains elements whose index values are the same as the
            exponent value that 10 is suposed to raised to.

    EXAMPLE:
        a = 1334.976
        resdiue(a)
        >>> [4,3,3,1]


    ATTENTION!:
        This function does not take into account the decimals of the input. (this was intentional yhough)

    THIS FUNCTION IS VERY BAD, IT HAS SOME MAJOR FLAWS AND BUGS THAT
    MAKE THE ANSWER RESULT WEIRD. BUT IT WORKS FOR SOME NUMBERS.

    """
    K = np.array([10**(-j) for j in range(40)])
    L = np.array([10**(j) for j in range(40)])

    power_index  = np.where(a > L)
    power_value = len(list(power_index[0]))
    power_decomp =[]
    absa = np.abs(a)
    for i in range(1,power_value+1):
        for j in range(1000):
            absa -= 10**(power_value-i)
            if absa < 0:
                power_decomp.append(j)
                absa = absa + 10**(power_value-i)
                break
            if absa == 0:
                power_decomp.append(0)
            if  1<= np.abs(absa) <= 2:
                power_decomp.append(1)
    l = len(power_decomp)
    pd = power_decomp
    val = 0
    for i in range(l):
        val += pd[i]*10**(i)



    return power_decomp[::-1], l, power_value,val

print(power_decomposition(3301601.8765))













"""
L = np.array([10**(j) for j in range(40)])
print(np.where(1000>L))
print(L)

"""





"""
    for i in range(n):
        res = a*10**(-i)
        ares = np.abs(res)
        exponent_index = list(np.where(ares > L))
        exponent_value = L[exponent_index]

        if ares < 10:
            break
    return res
print(residue(5.89))



def rounding_(a, decimal_point, n=100):
    decimal_strage =[]
    decimal_decomp =[]
    for i in range(n):
        res = a*10**(-i)
        if np.abs(res) < 10:
            break

    if res > 0:
        while res > 0:
            res -= 1
        decimal_storage.append(res)

    elif  res < 0:
        while res < 0:
            res += 1
        decimal_storage.append(res)

    mult_res = decimal_point*res
    k = decimal_point

    for i in range(k-1):
        for j in range(11):
            mult_res -= 10**(k-i)
            if mult_res < 0:
                decimal_decomp.append(j-1)
    return decimal_decomp, decimal_storage

print(rounding_(5.843, 3))
"""











def add(x,n):
   for i in range(n):
       x += y
   return x
print(add(0,4))









def fastlogapprox_amer(x,n, fast_diff_mat="off"):

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

X = np.linspace(2,10,10)
X = list(X)
V=[[fastlogapprox(x,s) for x in X ] for s in range(1,10)]
print(V)

for j in range(10):
    for x in X:
        u_new = fastlogapprox_amer(x,j)
        u.append(u_new)







def log_approx(n,x,difference = "off"):
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

    if difference == "off":
        return (x-1)/ a[-1]
    elif difference == "on":
        return np.abs((x-1)/a[-1] -np.log(x))
    elif difference == "fast":
        return  (x-1)/d[-1,-1]







print(log_approx(5,30))

X = np.linspace(2,100,100)
plt.figure()
for i in range(10):
    plt.plot(X, log_approx(i, X, difference ="fast"), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the error")






def fastlogapprox_hector(x, n):
    a = np.zeros(n)
    g = np.zeros(n)
    d = np.zeros((n,n))
    a[0] = 0.5*(1+x)
    g[0] = np.sqrt(x)
    for i in range(1,n):
        a[i] = 0.5*(a[i-1]+g[i-1])
        g[i] = np.sqrt(a[i]*g[i-1])
    d[0,:] = a
    for i in range(0,n):
        d[i,i] = (a[i] - 4**(-i)*a[i-1])/(1-4**(-i))
    return (x-1)/d[-1,-1]



def fastlogapprox_alex(x, n):
    a0 = 0.5*(1+x)
    g0 = np.sqrt(x)
    a = [a0]
    g = [g0]
    d = np.zeros((n,n))
    for i in range(n):
        if i<n:
            a.append(0.5*(a[i]+g[i]))
            g.append(np.sqrt(a[i+1]*g[i]))
            d[0,i] = a[i]
        elif i==n:
            break
    for k in range(1, n):
        for i in range(1, n):
            d[k,i] = ((d[k-1,i] - 4**(-k)*d[k-1, i-1])/(1-4**(-k)))
        return (x-1)/d[k,i]



def fastlogapprox_amer(x,n, fast_diff_mat="off"):

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
X = np.linspace(2,100,100)
u =[]
for i in range(10):
    u_new = fastlogapprox_amer(5,i)
for k in u:
    plt.plot(u[i],label=f"{i} iterations")
    plt.xlabel = ("x-value")
    plt.ylabel = ("y-value")
    plt.title("ln(x) and approx_ln(x)")
plt.show()


print(np.abs(fastlogapprox_hector(5, 3)))
print("-------------------------------------------------------------")
print(np.abs((fastlogapprox_alex(5, 3))))
print("-------------------------------------------------------------")
print(np.abs(fastlogapprox_amer(5, 3)))


X = np.linspace(2,100,100)
plt.figure()
for i in range(3):
    plt.plot(X, fastlogapprox_amer(X,i), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the function for different n")






print(np.abs(fastlogapprox_hector(5, 3)-np.log(5)))
print("-------------------------------------------------------------")
print(np.abs((fastlogapprox_alex(5, 3)-np.log(5))))
print("-------------------------------------------------------------")
print(np.abs(fastlogapprox_amer(5, 3)-np.log(5)))







def fastlogapprox_amer(n, x, fast_diff_mat="off"):

    a,g = [(1+x)/2], [np.sqrt(x)]

    for i in range(n):
        a_new = (a[-1] + g[-1])/2
        a.append(a_new)
        g_new = np.sqrt(g[-1]*a[-1])
        g.append(g_new)

    a = np.array(a)
    d = np.zeros((n,n))
    d[:1] = a
    for j in range(1,n):
        for i in range(j,n):
            d[j][i] = (d[j-1,i] -(d[j-1,i-1]*4**(-j)))/((1-4**(-j)))

    if fast_diff_mat == "on":
        d[-1:][0] = np.abs((x-1)/d[-1,-1] - (x-1)/a[-1])
        return d
    elif fast_diff_mat=="off":
        return d,(x-1)/d[-1,-1]




X = np.linspace(2,100,100)
plt.figure()
for i in range(1):
    plt.plot(X, fastlogapprox(X,i), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the function for different n")




# d[j-1:j, i-1:i] = ((d[j-1,i] - 4**(-j)*d[j-1,i-1]))/(1-4**(-j))








K = np.array([[1,2,8],
              [3,55,6],
              [7,8,89]])










x = np.linspace(2,5,100)


"""
HOMEWORK 1 lOG
"""
def log_approx(x, n, difference = "off", fast="on"):
    a,g ,f,t= [(1+x)/2], [np.sqrt(x)], [1], [1]
    d = np.zeros((n,n))

    for i in range(1,n):
        f_new = 4**(-i)
        t_new = 1-f_new
        f.append(f_new)
        t.append(t_new)

    for i in range(1,n):
        a_new = (a[-1] + g[-1])/2
        g_new = np.sqrt(g[-1]*a[-1])
        a.append(a_new)
        g.append(g_new)

    d[0] = a

    for j in range(1,n):
        for i in range(j+1, n):
            new = ((d[j,i] -f[j])**d[j,i-1])/t[j]

    if fast =="off":
        if difference == "off":
            return (x-1)/ a[-1]
        elif difference == "on":
            return (x-1)/a[-1] -np.log(x)

    elif fast =="on":
        if difference == "off":
            return (x-1)/ d[n+1,n+1], d
        elif difference == "on":
            return (x-1)/d[n,n] -np.log(x)
print(log_approx(5,10))








X = np.linspace(2,100,100)
plt.figure()
for i in range(5):
    plt.plot(X, log_approx(X,i, diff = "off"), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the function log_approx")





def s(x, *args):
    for s in args:
        x*=s
    return x

def functional(f,x,n, func_ =None):
    args = [k for k in range(n+1)]

    if func_ == "factorial":
        x = 1
        return s(x,*args)
    elif not func_:
        return s(x,*args)






"""
LECTURE 5

Unpacking arguments:
    positiona arhumets:
        - def f(a,b,c,d):
            return "some return"
        - what would happen if L =[1,2,3,4]
          and we wanted to printout f(L)?, well we would get an error.
          we should instead unpack the list using f(*L)
        - now, if the dictionary d={"b":1, "d":4, "c":3, "k":8}
          and do the following: f(d), we'd get an error, so we instead swrite
          f(**d)
    Passing (tunneling) argumnts:
       - def uter(f,x, *args, **keywords):
           return f(x,*args)

       - def f(x,g):
           return lambda x: g(x,g)
       - def newton(f, x0, tol, *args):
           ...
           f(x, *args)
    throw away varaiable:
        def my_function(x):
            return 1,2,3,4,5,6
        l = my_function(20)
        a = l[0]
        a,b,c,d = my_function
        a,b,_,d = my_functioon:
            - c is a "throw away variabl" and is not assigned any value













"""











def sign_separation_list_comprehension(A, content ="asnzp", neg =[], zero=[], pos=[], string=[]):
    """
    This function takes in a list as a positional arhument and returns a list
    containing some amount of lists (depending on the value of the deafult a
    rgument "content)

    ARGUMENT
    ----------
        A (list)
        content (str) can either be "asnzp", "snzp" or "nzp"


    RETURN
    ------
    (list)
    """
    neg    =[s for s in A if  s<0]
    zero   =[s for s in A if s== 0 ]
    pos    =[s for s in A if s>0 ]

    return [neg, zero, pos]

print(sign_separation_list_comprehension([1,2,3,4,5,6,-1.,0,0,-2,-3,-4,-5.,-6.]))












help(quad)

help(dblquad)

"""
HOMEWORK 1 lOG
"""
def log_approx(x, n, difference = "off", fast="on"):
    a,g ,f,t= [(1+x)/2], [np.sqrt(x)], [1], [1]
    d = np.zeros((n,n))

    for i in range(1,n):
        f_new = 4**(-i)
        t_new = 1-f_new
        f.append(f_new)
        t.append(t_new)

    for i in range(1,n):
        a_new = (a[-1] + g[-1])/2
        g_new = np.sqrt(g[-1]*a[-1])
        a.append(a_new)
        g.append(g_new)

    d[0] = a

    for j in range(1,n):
        for i in range(j+1, n):
            d[j][i] = ((d[j,i] -f[j])**d[j,i-1])/t[j]

    if fast =="off":
        if difference == "off":
            return (x-1)/ a[-1]
        elif difference == "on":
            return (x-1)/a[-1] -np.log(x)

    elif fast =="on":
        if difference == "off":
            return (x-1)/ d[n+1,n+1], d
        elif difference == "on":
            return (x-1)/d[n,n] -np.log(x)
print(log_approx(5,10))








X = np.linspace(2,100,100)
plt.figure()
for i in range(5):
    plt.plot(X, log_approx(X,i, diff = "off"), label = f"{i} iterations")

plt.legend(loc = "best")
plt.title("plotting the function log_approx")










"""
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
X = np.linspace(-500, 500, 1000)
plt.figure()
for k in range(5):
    plt.plot(X, poly(X, u[k]))
plt.title("noice plots")

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

X = np.linspace(-500, 500, 1000)
plt.figure()
plt.plot(X, PLY(X, c))
plt.title("noice plots")
"""


"""

"""
"""
HOMEWORK 1 lOG
"""
"""
def log_approx(x, n, diff = "on"):
    a,g = [(1+x)/2], [np.sqrt(x)]

    for i in range(n):
        a_new = (a[-1] + g[-1])/2
        a.append(a_new)
        g_new = np.sqrt(g[-1]*a[-1])
        g.append(g_new)
    if diff == "off":
        return (x-1)/ a[-1]
    elif diff == "on":
        return (x-1)/a[-1] -np.log(x)

x = np.linspace(2,5,100)

for i in range(5):
    plt.plot(X, log_approx(X,i))
plt.title("noice plots")


"""



































