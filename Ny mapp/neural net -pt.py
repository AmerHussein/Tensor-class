# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 15:50:58 2019

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

def levi_civita(rank, dim, detG, contravariant=None, covariant=None):
    if not contravariant and not covariant:
        raise Exception("you must choose if the levi-civita tensor either is co- or contravariant")
    else:
        def eps_gen(*args):
            def E(*args):
                if len(args)==2:
                    if args[0] > args[1]:
                        return False
                    elif args[0] < args[1]:
                        return True
                    else:
                        return None
                MIN = min(list(args))
                MAX = max(list(args))
                args = list(args)
                if MIN != 0:
                    comp = tuple([i for i in range(MIN, MAX+2)])
                    boool = [comp[i] in args for i in range(len(comp))]
                else:
                    comp = [i for i in range(MIN, MAX+1)]
                    boool = [comp[i] in args for i in range(len(comp))]

                if False in boool:
                    return None
                else:
                    args = list(args)
                    def LOL(*list_):
                        l = len(list_)
                        def lol(push_forward, *list_ , inv= False):
                            i = push_forward
                            E = list(range(l))
                            if inv == False:
                                L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
                            else:
                                L =  E[-(len(E) -(i+1)):] + E[:i+1]
                            return tuple(L)

                        e = [lol(i, *list_) for i in range(-1,l-1)]
                        r = [lol(i, *list_, inv=True) for i in range(-1,l-1)]
                        return e,r
                    Q = LOL(*args)
                    args = tuple(args)
                    if args in Q[1]:
                        return True
                    elif args in Q[0]:
                        return False
            a = E(*args)
            if contravariant:
                b = 1/np.sqrt(detG)
            elif covariant:
                b = np.sqrt(detG)
            if a==True:
                return b
            elif a ==False:
                return -b
            elif a ==None:
                return 0
        if rank == 2:
            eps = np.array([[eps_gen(*(i,j)) for i in range(dim)] for j in range(dim)])
        elif rank == 3:
            eps = np.array([[[eps_gen(*(i,j,k)) for i in range(dim)] for j in range(dim)] for k in range(dim)])
        elif rank ==4:
            eps = np.array([[[[eps_gen(*(i,j,k,u)) for i in range(dim)] for j in range(dim)] for k in range(dim)] for u in range(dim)])

        return eps

e = levi_civita(3,3, 1., covariant=True)
print(e)

def E(*args):
    if len(args)==2:
        if args[0] > args[1]:
            return False
        elif args[0] < args[1]:
            return True
        else:
            return None

    MIN = min(list(args))
    MAX = max(list(args))
    args = list(args)

    if MIN != 0:
        comp = tuple([i for i in range(MIN, MAX+2)])
        boool = [comp[i] in args for i in range(len(comp))]
    else:
        comp = [i for i in range(MIN, MAX+1)]
        boool = [comp[i] in args for i in range(len(comp))]

    if False in boool:
        return None
    else:
        args = list(args)
        def LOL(*list_):
            l = len(list_)

            def lol(push_forward, *list_ , inv= False):
                i = push_forward
                E = list(range(l))
                if inv == False:
                    L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
                else:
                    L =  E[-(len(E) -(i+1)):] + E[:i+1]
                return tuple(L)

            e = [lol(i, *list_) for i in range(-1,l-1)]
            r = [lol(i, *list_, inv=True) for i in range(-1,l-1)]
            return e,r

        Q = LOL(*args)
        args = tuple(args)
        if args in Q[1]:
            return True
        elif args in Q[0]:
            return False


R = E(*(0,1,2,3,4,5,6,7,7))

print(R)


def p(*args):
    if len(args) == 2:
        if args[0] > args[1]:
            return False
        elif args[0] < args[1]:
            return True
        else:
            return None

    MIN = min(list(args))
    MAX = max(list(args))
    args = list(args)

    if MIN != 0:
        comp = tuple([i for i in range(MIN, MAX+2)])
        boool = [comp[i] in args for i in range(len(comp))]
    else:
        comp = [i for i in range(MIN, MAX+1)]
        boool = [comp[i] in args for i in range(len(comp))]

    if False in boool:
        return None
    else:
        args = list(args)
        def LOL(*list_):
            l = len(list_)

            def lol(push_forward, *list_ , inv= False):
                i = push_forward
                E = list(range(l))
                if inv == False:
                    L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
                else:
                    L =  E[-(len(E) -(i+1)):] + E[:i+1]
                return tuple(L)

            e = [lol(i, *list_) for i in range(-1,l-1)]
            r = [lol(i, *list_, inv=True) for i in range(-1,l-1)]
            return e,r

        Q = LOL(*args)
        args = tuple(args)
        if args in Q[1]:
            return True
        elif args in Q[0]:
            return False



R = E(*(1,2,3,4,5))

print(R)






def index_perm(*args):
    """
    This function checks if a given tuple or list or numbers (integers)
    is a positive (returning True), negative (t´retuning False) or None of them
    retuning None.
    """
    l = len(args)

    args = list(args)
    MIN = min(args)
    MAX = max(args)
    comp = [i for i in range(MIN,  MAX+2)]

    boool = [comp[i] in args for i in range(len(comp))]
    if False in boool:
        return None
    else:
        def lol(push_forward, *args, inv= False):
            i = push_forward
            E = list(args)
            if inv == False:
                L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
            else:
                L =  E[-(len(E) -(i+1)):] + E[:i+1]
            return tuple(L)

        e = [lol(i, *args) for i in range(-1,l-1)]              #inverse
        r = [lol(i, *args, inv=True) for i in range(-1,l-1)]    #normal

        if args in e:
            return False
        elif args in r:
            return True
print(index_perm(*(1,2,3)))


def positive_perm(*args):
    """
    This function checks if a given tuple or list or numbers (integers)
    is a positive (returning True), negative (t´retuning False) or None of them
    retuning None.
    """
    l = len(args)

    args = list(args)
    MIN = min(args)
    MAX = max(args)
    """
    if MIN != 0:
        comp = [i for i in range(MIN,  MAX+1)]
        boool = [comp[i] in args for i in range(len(comp))]
    else:
        comp = [i for i in range(MAX+2)]
        boool = [comp[i] in args for i in range(len(comp))]

    if False in boool:
        return None
    else:
    """
    def lol(push_forward, *args, inv= False):
        i = push_forward
        E = args
        if inv == False:
            L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
        else:
            L =  E[-(len(E) -(i+1)):] + E[:i+1]
        return tuple(L)

    e = [lol(i, *args) for i in range(-1,l-1)]              #inverse
    r = [lol(i, *args, inv=True) for i in range(-1,l-1)]    #normal

    if args in e:
        return False
    elif args in r:
        return True


print(positive_perm(*(0,1,2)))

def levi_civita(rank, dim, detG, contravariant=None, covariant=None):
    if not contravariant and not covariant:
        raise Exception("you must choose if the levi-civita tensor is co- or contravariant")
    else:
        def eps_gen(*args):
            def E(*args):
                MIN = min(list(args))
                MAX = max(list(args))
                args = list(args)



                if MIN != 0:
                    comp = tuple([i for i in range(MIN, MAX+2)])
                    boool = [comp[i] in args for i in range(len(comp))]
                else:
                    comp = [i for i in range(MIN, MAX+1)]
                    boool = [comp[i] in args for i in range(len(comp))]

                if False in boool:
                    return None
                else:
                    args = list(args)
                    def LOL(*list_):
                        l = len(list_)
                        def lol(push_forward, *list_ , inv= False):
                            i = push_forward
                            E = list(range(l))
                            if inv == False:
                                L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
                            else:
                                L =  E[-(len(E) -(i+1)):] + E[:i+1]
                            return tuple(L)

                        e = [lol(i, *list_) for i in range(-1,l-1)]
                        r = [lol(i, *list_, inv=True) for i in range(-1,l-1)]
                        return e,r

                    Q = LOL(*args)
                    args = tuple(args)
                    if args in Q[1]:
                        return True
                    elif args in Q[0]:
                        return False


                a = index_perm(*args)
                if a==True:
                    return 1.
                elif a ==False:
                    return -1.
                else:
                    return 0

        if rank == 2:
            eps = np.array([[eps_gen(*(i,j)) for i in range(dim)] for j in range(dim)])
        elif rank == 3:
            eps = np.array([[[eps_gen(*(i,j,k)) for i in range(dim)] for j in range(dim)] for k in range(dim)])
        elif rank ==4:
            eps = np.array([[[[eps_gen(*(i,j,k,u)) for i in range(dim)] for j in range(dim)] for k in range(dim)] for u in range(dim)])

        if contravariant:
            return eps/np.sqrt(detG)
        elif covariant:
            return eps*np.sqrt(detG)
print(levi_civita(2, 4, 1., contravariant=True))





def lol(push_forward):
    i = push_forward
    E = [1,2,3,4,5]
    L =  E[-(len(E) -(i+1)):] + E[:i+1]
    return L
#for i in range(-1,4):
 #   print(lol(i))


def positive_perm(*args):
    """
    This function checks if a given tuple or list or numbers (integers)
    is a positive (returning True), negative (t´retuning False) or None of them
    retuning None.
    """
    l = len(args)

    args = list(args)
    MIN = min(args)
    MAX = max(args)
    comp = [i for i in range(MIN,  MAX+1)]

    boool = [comp[i] in args for i in range(len(comp))]
    if False in boool:
        return None
    else:
        def lol(push_forward, *args, inv= False):
            i = push_forward
            E = args
            if inv == False:
                L =  E[:i+1][::-1] + E[-(len(E) -(i+1)):][::-1]
            else:
                L =  E[-(len(E) -(i+1)):] + E[:i+1]
            return L

        e = [lol(i, *args) for i in range(-1,l-1)]              #inverse
        r = [lol(i, *args, inv=True) for i in range(-1,l-1)]    #normal

        if args in e:
            return False
        elif args in r:
            return True





print(positive_perm(*(1,2,3, 4)))



















































