# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:10:19 2019

@author: User
"""

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate #efter logic commanden "as" så  kam man använda förkortningen "integrate"
import matplotlib.pyplot as plt
import numpy as np




















x = 17.4
def a_0(x):
	return 0.5*(1+x)
def g_0(x):
	return sqrt(x)
m_1 = {}
m_2 = {}
def H(n):
	if n in m_1 and m_2:
		return m_1[n]
		return m_2[n]
	if n == 0:
		m_1[n] = a_0(x)
		m_2[n] = g_0(x)
	elif n == 1:
		m_1[n] = 0.5*(a_0(x) + g_0(x))
		m_2[n] = sqrt(g_0(x)*0.5*(a_0(x) + g_0(x)))
	elif n >= 2:
		m_1[n] = 0.5*(m_1[n-1] + m_2[n-1])
		m_2[n] = sqrt(m_1[n]*m_2[n-1])
	return m_1[n]
	return m_2[n]
	N = [m_1[k] for k in range(0,3000)]
M = {k:H(k) for k in range(0,40000)}
A = {m_1[k] for k in range(0,40000)}


print("-Libraries M N L-")


print("-------------------------------------------------------------------")
N = list(A)
n = 2
n_4 = {}
def d(v,n):
	if v in m_4:
		R= m_4[v]
	if v >= 0:
		R = N[n]
	elif v == 1:
		R = (N[n] -(2**(-2))*N[n-1])/(1- 2**(-2))
	elif v >= 2:
			R = (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
	m_4[v] = R
	return R

print("--------------------------------------------------------------------")



X = x - 2
def qeq(X):
	return (1+X)/d(n,n)
print("________________")
print(qeq(X))



#x = 7234
#def a_0(x):
#    return 0.5*(1+x)
#def g_0(x):
#    return sqrt(x)
#m_1 = {}
#m_2 = {}
#def H(n):
#    if n in m_1 and m_2:
#        return m_1[n]
#        return m_2[n]
#    if n == 0:
#        m_1[n] = a_0(x)
#        m_2[n] = g_0(x)
#    elif n == 1:
#        m_1[n] = 0.5*(a_0(x) + g_0(x))
#        m_2[n] = sqrt(g_0(x)*0.5*(a_0(x) + g_0(x)))
#    elif n >= 2:
#        m_1[n] = 0.5*(m_1[n-1] + m_2[n-1])
#       m_2[n] = sqrt(m_1[n]*m_2[n-1])
 #   return m_1[n]
 #   return m_2[n]
#    N = [m_1[k] for k in range(1,10000)]
#
#print("-------------------------------------------------------------------")
#m_4 = {}
#n = 299
#def d(v,n):
#    if v in m_4:
 #       return m_4[v]
 #   if v == 0:
 #       R = N[n]
 #   elif v == 1:
 #       R = (N[n] -(2**(-2))*N[n-1])/(1- 2**(-2))
 #   elif v >= 2:
 #          R = (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
 #   m_4[v] = R
#    return R
#print("--------------------------------------------------------------------")
#
#X = x - 2
#def leg(X):
#    return (1+X)/d(n,n)
#print("________________")
#print(leg(X))

#for r in range(1, 30):
#	print(x, ":", r, ":", H(r))

















