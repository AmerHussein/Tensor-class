# -*- coding: utf-8 -*-
"""
Created on Sun Apr 21 18:55:02 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate #efter logic commanden "as" så  kam man använda förkortningen "integrate"
import matplotlib.pyplot as plt
import numpy as np


x = 10
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
	N = [m_1[k] for k in range(0,30)]
print("-------------------------------------------------------------------")
m_4 = {}
n = 23
def d(v,n):
	if v in m_4:
		return m_4[v]
	if v == 0:
		R = N[n]
	elif v == 1:
		R = (N[n] -(2**(-2))*N[n-1])/(1- 2**(-2))
	elif v >= 2:
			R = (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
	m_4[v] = R
	return R
print("--------------------------------------------------------------------")

X = x - 2
def leg(X):
	return (1+X)/d(n,n)
print("________________")
print(leg(X))


	#for r in range(1, 10000):
	#print(x, ":", r, ":", H(r))

x = 8000000000000000000
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
    return m_1[n] and m_2[n]


Q = {r:H(r) for r in range(1, 40000)}

print("-Libraries M N L-")
M = [H(k) for k in range(0,40000)]
A = [m_1[k] for k in range(0,40000)]
#print("--------------------------------------------")
m_4 = {}
N = list(A)
n = 2950
def d(v,n):
    if v in m_4:
        return m_4[v]
    if v == 0:
        R = N[n]
    elif v == 1:
        R = (N[n] -(2**(-2))*N[n-1])/(1- 2**(-2))
    elif v >= 2:
        R = (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
    m_4[v] = R
    return R
#print("-------------------------------------------")

X = x-2
def leg(X):
    return (1+X)/d(n,n)
print("______________")
print(leg(X))



#x = np.linspace(1, 20)
#plt.plot(x, (leg(x)),color='blue', label='approx ln(x)' )
#plt.plot(x, np.log(x), color='red', label='ln(x)')
#plt.legend(loc='upper left')
#plt.show()
#plt.plot(x,np.abs(np.log(x) - leg(X)) )
#plt.show()

plt.plot(x,np.abs(np.log(x) - leg(x)) )
plt.show()



























