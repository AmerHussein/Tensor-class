# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 20:46:57 2019

@author: User
"""







def zero(f, a, b, margin=1.e-2, epsilon=1.e-3, max_iter_ = 30000):
    def rec(f, a,b, c= "linspace", div= 100,  max_iter=300, step=0.01, tol=1.e-3):
            z =[(a+b)/2, ((a+b)/2 +b)/2]
            if c == "usual":
                for k in range(1, max_iter):
                    if z[-1]-z[-2] > tol:
                        z.append((z[k-1]+b)/2)
                    else:
                        z.append(z[k-1] + step)
                        if z[-1] > b:
                            break
                return z
            elif c == "linspace":
                F = np.linspace(a,b, div)
                R = np.zeros_like(F)
                for k in range(len(R)):
                    R[k] = f(F[k])
                return R

    def mainframe(f, rec):
        R = rec(a,b)
        G = np.zeros_like(R)
        l = len(list(R))
        if f(a)>0>f(b):
            for k in range(l):
                G[k] = f(R[k])
                if G[k]==0 or np.abs(G[k]) <= margin:
                    break
                    return R[k]
                elif G[k]<0:
                    Z = np.zeros(l)
                    q = np.zeros(l)
                    Z[k] = R[k] - epsilon
                    q[k] = f(Z[k])
                    if q[k]==0 or np.abs(q[k]) <= margin:
                        break
                        return Z[k]
                    elif q<0:
                        Q = np.zeros(max_iter_)
                        for j in range(max_iter_):
                            Q[j] = f(Z[k]-j*epsilon)
                            if np.abs(Q[j])<= margin:
                                break
                        return Z[k]-j*epsilon
        elif f(b)>0>f(a):
            R = rec(a,b)
            G = np.zeros_like(R)
            l = len(list(R))
            if f(a)>0>f(b):
                for k in range(l):
                    G[k] = f(R[k])
                    if G[k]==0 or np.abs(G[k]) <= margin:
                        break
                        return R[k]
                    elif G[k]>0:
                        Z = np.zeros(l)
                        q = np.zeros(l)
                        Z[k] = R[k] - epsilon
                        q[k] = f(Z[k])
                        if q[k]==0 or np.abs(q[k]) <= margin:
                            break
                            return Z[k]
                        elif q>0:
                            Q = np.zeros(max_iter_)
                            for j in range(max_iter_):
                                Q[j] = f(Z[k]-j*epsilon)
                                if np.abs(Q[j]) <= margin:
                                    break
                            return Z[k]-j*epsilon

        elif f(a)==0 or np.abs(f(a))<= margin:
            return a
        elif f(b)==0 or np.abs(f(b))<= margin:
            return b



def p(x):
    return -x +6
print(zero(np.cos,0,np.pi))























def rec(a,b, c= "linspace", div= 100,  max_iter=30, step=0.01, tol=1.e-3):
        z =[(a+b)/2, ((a+b)/2 +b)/2]
        if c == "usual":
            for k in range(1, max_iter):
                if z[-1]-z[-2] > tol:
                    z.append((z[k-1]+b)/2)
                else:
                    z.append(z[k-1] + step)
                    if z[-1] > b:
                        break
            return z
        elif c == "linspace":
            R = np.linspace(a,b, div)
            return R



print(rec(0,12))


def zero(f, a, b, margin=1.e-4, epsilon=1.e-5, linspace_divider=100):

    def rec(a,b, c= "linspace", div= 100,  max_iter=30, step=0.01, tol=1.e-3):
        z =[(a+b)/2, ((a+b)/2 +b)/2]
        if c == "usual":
            for k in range(1, max_iter):
                if z[-1]-z[-2] > tol:
                    z.append((z[k-1]+b)/2)
                else:
                    z.append(z[k-1] + step)
                    if z[-1] > b:
                        break
            return z
        else:
            R = np.linspace(a,b, div)
            return R

    z = rec(a,b)
    l = len(z)



    def mainframe(f, z):
        for i in range(l):
            F = f(z[i])
            if np.abs(F) <= margin or F==0:
                break
            return z[i]

            if F < 0:
                Z = z[i] - epsilon
                q = f(Z)
                if q < 0:
                    for j in range(l-i),l:
                        Q = f(Z -j*epsilon)
                        if np.abs(Q) <= margin or Q ==0:
                            break
                    return Z-j*epsilon


    if f(a)>0>f(b):
        yield mainframe(f,z)
        if mainframe(f,z) is None:

            x = np.linspace(a,b, np.linspace_divider)
            X = len(x)
            for k in range(X):
                F = f(x[k])
                if F ==0 or np.abs(F)<= margin:
                    break
                    return x[k]
                elif F<0:
                    E = rec(x[k],b)
                    return mainframe(f, E)



    #conditional preparation in case the input gives none or an error.

    elif f(a)>0 and f(b)>0:
       x = np.linspace(a,b, np.linspace_divider)
       X = len(x)
       for k in range(X):
          F = f(x[k])
          if F ==0 or np.abs(F)<= margin:
              break
              return x[k]
          elif F<0:
              E = rec(x[k],b)
              return mainframe(f, E)

#    elif f(a)<0<f(b):
#        f =
#        return print(mainframe(f,z))


def p(x):
    return -x +6
print(zero(p,0,12))













































def zero(f, a,b , max_iter =15):

    z = (a+b)/2
    for k in range(max_iter):
        if  f(a)*f(b) > 0:
            #raise Exception("EVT does not hold")
            continue

        elif f(a)>0>f(b):
            a = (b +z)/2
            a = z

        elif f(b)>0>f(a):
            b = (a+z)/2
            b = z
    return z

def p(x):
    return x -6
print(zero(p,0,10))


def f(a,b,n):
    if n==0:
        return a
    elif n == 1:
        return (a+b)/2
    elif n ==2:
        return (a+b)/2**2
    elif n > 1:
        return f(a,b,n-1)

print(f(1,2,100))






def zero(f, a,b, max_iter =10):
    def recursive(a,b,n):
        if n ==1:
            return (a+b)/2
        elif n>1:
            return recursive(a,b, n-1)
    for k in range(1,max_iter):
        if f(b)>0>f(a):
            m = recursive(a,b, k)
            if f(m)<0:
                m = a
                b = b
            elif f(m)>0:
                a = a
                m = b
            elif f(m)==0:
                break
                return m
        #elif f(b)<0<f(a):

    return m


def p(x):
    return x -6
print(zero(p,0,10))





def div_n_conc(f,a,b, tol= 1.e-10, max_iterk =10, max_iterK= 10, k=0, K=0):
    z=[]
    z.append((a+b)/2)
    q = np.array([[a,z[-1]], [z[-1],b]])

    while K<max_iterK:
        if f(b)>0>f(a) and np.abs(f(b)) > np.abs(f(a)):
            while k<max_iterk:
                b = z[k]
                z[k+1] += (a+z[k])/2
                if not f(b)>0>f(a) and not np.abs(f(b)) > np.abs(f(a)):
                    break
                else:
                    k+=1
        else:
            k = 0

        if f(b)>0>f(a) and np.abs(f(b)) < np.abs(f(a)):
            while k<max_iterk:
                a = q[1,0]
                a += (a+b)/2
                if not f(b)>0>f(a) and not np.abs(f(b)) < np.abs(f(a)):
                    break
                else:
                    k+=1
        else:
            k = 0

        if f(a)>0>f(b) and np.abs(f(a)) > np.abs(f(b)):
            while k<max_iterk:
                a = q[1,0]
                a += (a+b)/2
                if not f(a)>0>f(b) and not np.abs(f(a)) > np.abs(f(b)):
                    break
                else:
                    k+=1
        else:
            k =0

        if f(a)>0>f(b) and np.abs(f(b)) > np.abs(f(a)):
            while k<max_iterk:
                b = q[0,1]
                b += (a+b)/2
                if not f(a)>0>f(b) and not np.abs(f(b)) > np.abs(f(a)):
                    break
                else:
                    k +=1
        K +=1
        return q

def p(x):
    return x -3
print(div_n_conc(p, 0,4))
















def div_n_conc(f,a,b,C = "short answer", tol= 1.e-10, max_iterk =5, max_iterK= 100, k=0, K=0):
    z = np.zeros(max_iterk)
    z[0] = (a+b)/2
    q = np.array([[a,z[0]], [z[0],b]])

    if C == "short answer":
#        while K < max_iterK and np.abs(q[1][0] -q[0][1]) < tol:
        while K < max_iterk and np.abs(q[1][1] -q[0][1]) > tol:
            while f(b)>0>f(a) and np.abs(f(b)) > np.abs(f(a)):
                z[k+1]  = (z[k] +b)/2
            break
            k +=1
            if f(b)>0>f(a) and np.abs(f(b)) < np.abs(f(a)):
                b = (q[1,0]+a)/2


            elif f(a)>0>f(b) and np.abs(f(a)) > np.abs(f(b)):
                b = (q[1,0]+a)/2


            elif k < max_iterk and f(a)>0>f(b) and np.abs(f(b)) > np.abs(f(a)):
                a = (q[1,0]+b)/2
                k += 1
        K += 1
        return q
    elif C == "long answer":
        return 2

def p(x):
    return x -3


print(div_n_conc(np.cos, 1, np.pi))




