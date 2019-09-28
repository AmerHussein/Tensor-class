# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 13:24:03 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np



def fastlogapprox(x, n):
    a = zeros(n)
    g = zeros(n)
    d = zeros((n,n))
    a[0] = 0.5*(1+x)
    g[0] = sqrt(x)
    for i in range(1,n):
        a[i] = 0.5*(a[i-1]+g[i-1])
        g[i] = sqrt(a[i]*g[i-1])
    d[0,:] = a
    for i in range(0,n):
        d[i,i] = (a[i] - 4**(-i)*a[i-1])/(1-4**(-i))
    return (x-1)/d[-1,-1]

for k in range(1, 5):
    x_vals = np.linspace(0, 100, 1000)
    y_vals = fastlogapprox(x_vals, k)
    plt.plot(x_vals, y_vals)
    plt.title ("shit")
    plt.xlabel('x')
    plt.ylabel('square root of x')
    plt.legend()
    plt.show()



q = np.abs(np.log(3) - fastlogapprox(3,2))
w = np.abs(np.log(3) - fastlogapprox(3,3))
r = np.abs(np.log(3) - fastlogapprox(3,4))
t = np.abs(np.log(3) - fastlogapprox(3,5))
u = 1e-19
z = 1e-05

plt.plot(3, q, color='blue', label='2 Iterations')
plt.plot(3, w, color='green', label='3 Iterations')
plt.plot(3, r, color='red', label='4 Iterations')
plt.plot(3, t, color='cyan', label='5 Iterations')
plt.legend(loc='upper left')
plt.title('Error behavior of the accelerated Carlsson Method for the log')
plt.xlabel('x')
plt.ylabel('error')
plt.axis([0, 20, u, z])
plt.show()













#def applog(a,x):
#    a0=(1+x)/2
#   b0=sqrt(x)
#    a=[]
#
#   for i in range(n):
#        a0=(a0+b0)/2
#        b0=sqrt((a0)*b0)
#        a.append(n)
#
#    return np.array((x-1)/a0)




#print("the log approximation is: " , ((applog(4,4))))
#print("log value----------------:" , (np.log(4)) )

#error=(abs( applog(4,4)- (np.log(4))   ) )
#print("the error is : ", error)


#x = np.linspace(10, 500, 100)
#plt.plot(x, (applog(4,x)),color='blue', label='approx ln(x)' )
#plt.plot(x, np.log(x), color='red', label='ln(x)')
#plt.legend(loc='upper left')
#plt.show()
#plt.plot(x,abs(np.log(x)-applog(4,x)) )
#plt.show()













print("----------------------------------------------------------------------")












#x = 4
#def a_0(x):
#	return 0.5*(1+x)
#def g_0(x):
#	return sqrt(x)
#m_1 = {}
#m_2 = {}
#def H(n):
#	if n in m_1 and m_2:
#		return m_1[n]
#		return m_2[n]
#	if n == 0:
#		m_1[n] = a_0(x)
#		m_2[n] = g_0(x)
#	elif n == 1:
#		m_1[n] = 0.5*(a_0(x) + g_0(x))
#		m_2[n] = sqrt(g_0(x)*0.5*(a_0(x) + g_0(x)))
#	elif n >= 2:
#		m_1[n] = 0.5*(m_1[n-1] + m_2[n-1])
#		m_2[n] = sqrt(m_1[n]*m_2[n-1])
#	return m_1[n] and m_2[n]


#Q = {r:H(r) for r in range(1, 40000)}
#M = [H(k) for k in range(0,40000)]
#A = [m_1[k] for k in range(0,40000)]
#print("--------------------------------------------")
#m_4 = {}
#N = list(A)
#n = 3
#X = x-2
#def leg(x,n):
#	def d(v,n):
#		if v in m_4:
#			return m_4[v]
#		if v == 0:
#			R = A[n]
#		elif v == 1:
#			R = (A[n] -(2**(-2))*A[n-1])/(1- 2**(-2))
#		elif v >= 2:
#			R = (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
#		m_4[v] = R
#		return R
#	return (1+x)/d(n,n)
#print(leg(X,n))
#
#x_v =  linspace(3, 4003, 100000)
#n_q = range(1, 10)
#y_1 = leg(2,n)
#y_2 = abs(leg(2,n_q)- log(2))
##plot(x_v ,y_1, label='leg(x,n)')
#plot(x_v ,y_2, label='error')
#title ('plopp')
#xlabel('x')
#ylabel('leg(x)')
#legend()
#show()



#print("-------------------------------------------")

#X = x-2
#def leg(X):
#	return (1+X)/d(n,n)
#print("______________")
#print(leg(X))















