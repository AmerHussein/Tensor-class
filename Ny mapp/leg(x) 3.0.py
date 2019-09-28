# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 14:12:46 2019

@author: User
"""


from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np







x = 1.6
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
n = 287
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


def leg(X):
	return (1+X)/d(n,n)
print("______________")
print(leg(X))

d = []
def b(n):
	for n in range(0,100):
		d.append(A[n])





x_v =  linspace(3, 4003, 100000)
y_1 = leg(x_v)
y_2 = abs(leg(x_v)- log(x_v))
plot(x_v ,y_1, label='leg(x)')
plot(x_v ,y_2, label='error')
title ('plopp')
xlabel('x')
ylabel('leg(x)')
legend()
show()






def logapprox(x, n):
    a = zeros(n)
    g = zeros(n)
    a[0] = 0.5*(1+x)
    g[0] = sqrt(x)
    for i in range(1,n):
        a[i] = 0.5*(a[i-1]+g[i-1])
        g[i] = sqrt(a[i]*g[i-1])
    return (x-1)/a[-1]
print(logapprox(80, 310))






































