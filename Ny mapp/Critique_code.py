# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 10:34:39 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import scipy.linalg as sl

"""
---------------------------LOOK BELOW THE FUNCTION---------------------
"""


"¯\_(ツ)_/¯"


def Beta_bool(B,n, max_iter =30, margin = 1.e-9, x0=1):
    """
    denna funktion genomför n rekusrion under en for-loop
    ON INPUT:
        B (float)
    """
    c= [0.2  ,5]
    x=[x0]
    for i in range(n):
        x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
        x.append(x_new)
        delta = np.abs(x[-1]- x[-2])
        if i < max_iter:
            if delta < margin:
                return True, x

        elif i == max_iter:
            if delta > margin:
                return False, x

print(Beta_bool(0.05, 100))


def Beta(B,n, p = "x", margin = 1.e-9, x0=1):
    """
    this function executes a recursion with the aid of a for-loop
    ON INPUT:
        *B (float):
            it is very important that we keep the star.
    ON OUTPUT:
        x (list)
        or
        x[-1] (float)
        or
        conv: {x[-1]} (string)

    """
    c,x = [0.2  ,5], [x0]
    for i in range(n):
        x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
        x.append(x_new)
        delta = np.abs(x[-1]- x[-2])
        if delta < margin:
            break

    if p == "x_lastval_sring":
        return f"""conv: {x[-1]}"""

    elif p == "x_lastval":
        return x[-1]

    elif p == "x":
        return x
print(Beta(0.5, 100))



def sign_separation(A, neg=[], zero=[], pos=[]):
    for s in A:
        if s > 0:
            pos.append(s)
        elif s < 0:
            neg.append(s)
        else:
            zero.append(s)

    if len(zero) == 0:
        return [neg, pos]
    else:
        return [neg, zero, pos]

b = [-0.5, 0.5, -0.25, 0, 0.25]
print(sign_separation(b))





def seq_conv(B,n, *var, p = "x_lastval",convergence_bool_swich = "on", sign_separation_swich ="off", bool_convert_swich = "on", margin = 1.e-9, x0 =1):

    """
    This function takes a list; var, and iterates the function Beta over its elements.
    You can shoose if you want to activate its different subfunctions, such as bool_convert,
    sign_separation or of you want to analyze the convergence of the recursion within some maximum
    iteration count.


    ON INPUT:
        n (int):
            amount of itertions
        var (list)
        p (string)

    ON OUTPUT:
        bb (list)
    """

    bb = []
    for s in var:
        s = Beta(s,n)
        bb.append(s)

    if bool_convert_swich == "on":
        V = []
        for s in var:
            v = Beta_bool(s,n)
            V.append(v)
        return V


    if sign_separation_swich == "on":
            return sign_separation(bb)
    elif sign_separation_swich == "off":
        return bb

    if convergence_bool_swich == "on":
        return Beta_bool(B, n)


b = [-0.5, 0.5, -0.25, 0, 0.25]
print(seq_conv(100, *b))






"""
--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
----The function at the top works fine i guess (not sure though)----------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
--------The function at the bottom does not work, i think-----------------
--------------------------------------------------------------------------
--------------------------------------------------------------------------
"""






def seq_conv(n, *var, p = "x_lastval", sign_separation_swich ="off", iteration_regulation = "on", margin = 1.e-9, x0 =1):
    """
    Denna funtion tar en lista *var, och genomför funcktionen Beta över dess element

    ON INPUT:
        n (int):
            amount of itertions
        var (list)
        p (string)

    ON OUTPUT:
        bb (list)

    """
    x,c,b = [x0], [0.2  ,5], []




    if iteration_regulation == "on":
        def Beta_moderator(B, max_iter = 30):
            """
            denna funktion genomför n rekusrion under en for-loop
            ON INPUT:
                B (float)
            """
            for i in range(max_iter+1):
                x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
                x.append(x_new)
                if np.abs(x[-1]- x[-2])< margin:
                    break
                    return True
                else:
                    return False

    elif   iteration_regulation == "off":
        def Beta(B):
            """
            denna funktion genomför n rekusrion under en for-loop
            ON INPUT:
                B (float)
            """
            for i in range(n):
                x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
                x.append(x_new)
                if np.abs(x[-1]- x[-2])< margin:
                    break

            if p == "x_lastval_sring":
                return f"""conv: {x[-1]}"""

            elif p == "x_lastval":
                return x[-1]

            elif p == "x":
                return x

    """
    here we apply Beta to our input var:
    """
    for s in var:
        s = Beta(s)
        b.append(s)


    if sign_separation_swich == "yes":
        """
        this functionseparates the elemts of list, A, into a list containing three lists,
        with negative, zero elements and positive elemts respectivelly
        """
        def sign_separation(A, neg=[], zero=[], pos=[]):
            for s in A:
                if s > 0:
                    pos.append(s)
                elif s < 0:
                    neg.append(s)
                else:
                    zero.append(s)

                if len(zero) == 0:
                    return [neg, pos]
                else:
                    return [neg, zero, pos]
        return sign_separation(b)

    elif sign_separation_swich == "no":
        return b





B = [-0.5, 0.5, -0.25, 0.25]
print(seq_conv(100, *B))































