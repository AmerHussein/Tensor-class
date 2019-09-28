# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 14:35:54 2019

@author: User
"""
from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate #efter logic commanden "as" så  kam man använda förkortningen "integrate"
import matplotlib.pyplot as plt
import numpy as np
from decimal import *
getcontext().prec = 25
#i = 2
#j = 1
#k = 3

#A = ['','smootie', 'capochino', 'black coffe']
#B = ['', 'milk', 'sugar', 'cream']
#C = ['', 'muffins', 'kex', 'cake']
#A.append(B)
#A.append(C)
#A =[['', 'smootie', 'capochino', 'black coffe'], ['', 'milk', 'sugar', 'cream'], ['', 'muffins', 'kex', 'cake']]
#i = A[0][i]
#j = A[1][j]
#k = A[2][k]
#print('your order is', i, 'with', j, 'and', k)

#exec(open('kemans kod.py').read(0))

#-------SIMPLE MODULES-COLLECTING FUNCTIONS------
# often oen collects fucntions in a script. this creates a
# module with additional python functionality. To demosntrate this, we create
# a module by collecting functions n a single file. for example kemans kod.py

def f(x):
	return 2*x + 1
def g(x): 
	return x**2
def h(x):
	return g(f(x))
for a in range(1, 5):
	if h(a) > 1:
		exec(open('FIZZBUZZ.py').read(0))
	else:
		print('none') 

# these functions can now be used by any any external script or directly in the 
# Ipyhton environment som tex i din kod för fizzbuz där du definierade ett antal funktioner och kune sedan använda dessa funktioner som då bildade en modul för resten av din kod. 
# funktioner i modulen kan bero på varandra
# OM MAN SKRIVER IN EN MODUL I EN HELT ENSKILD FIL OCH VET BETENDET AV MODULEN SÅ KAN VI KALLA IN DEN 
# I EN SCRIPT I EN ANNAN FIL. VI TESTAR
		#fuck it
		
#alternativlly, the modules can be imported by the command "import". It creates 
#  a named namspace. The command "from" puts the functions in the modules into 
# a general namespace.

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
for r in range(1, 30):
	print(x, ":", r, ":", H(r))
print("-Libraries M N L-")
M = {k:H(k) for k in range(0,11)}
print(M)

#print("--")

N = [m_1[k] for k in range(0,22)]
print(N)
for k in range(4):
	print(N[k])
n = 20
def d(v,n):
	if v == 0:
		return N[n]
	elif v == 1:
		return (N[n] -(2**(-2))*N[n-1])/(1- 2**(-2))        
	elif v >= 2:
			return (d(v-1, n) -(2**(-2*v))*d(v-1, n-1))/(1-2**(-2*v))
for k in range(1, 20):
	print(d(k, n))

print("----")

P = [d(k, n) for k in range(0,21)]


X = x-2
def leg(X):
	return (1+X)/d(n,n)
print("______________")
print(leg(X))
	


	
	
	


	
	










































