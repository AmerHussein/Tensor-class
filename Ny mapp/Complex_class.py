# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 22:42:29 2019

@author: User
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May  8 14:50:35 2019

@author: User
"""
"""
Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap,
 CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r,
 Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r,
 PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r,
 RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1,
 Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r,
 YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone,
 bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r,
 cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat,
 gist_heat_r,
 gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg,
 gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r,
 inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean,
 ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic,
 seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b,
 tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
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

x, y = symbols('x y', real=True)
f = y*(x/(1+x)) +1 +x
F1 = 3*x**3-3*x*y**2 -1
F2 = 3*x**2*y -y**3

a,b= diff(F1,x), diff(F1,y)
c,d= diff(F2, x), diff(F2, y)

print(a)
print(b)
print(c)
print(d)


d_x = diff(f, x)
d_y = diff(f, y)

print(d_x)
print(d_y)


def g(x):
    return x**2


f = lambdify((x,y),g(x), "numpy")
print(f(1, 2))



def motor_integration(a,b,c,d,N, F, J, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT=False, zero=[], plotting_method="mine", coloroption="plasma"):
    def tensor_motor(a,b,c,d,N):
        x = np.linspace(a,b,N)
        y = np.linspace(c,d,N)
        L = []
        X,Y = np.meshgrid(x,y)

        for i in range(N):
            L.append([])

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                L[i].append((X[i,j], Y[i,j]))
        L = np.array(L)
        return  L

    def Newt_motor(F, J, v0, max_iter=max_iter, margin=1.e-4, step=step, h=1.e-9, finite=False, once=False):
        v0 = np.array(v0)
        xvals = [v0]
        V = xvals[-1]
        V_nrm = np.linalg.norm(xvals[-1])
        F_nrm = np.linalg.norm(F(xvals[-1]))
        x_max = 1/margin
        F_max = 1/margin
        k = 0

        def solving_motor(v):
            def finite_difference(v, h):
                J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                              [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                return J
            if finite == False:
                if once==False:
                    delta = np.linalg.solve(J(v), -F(v))
                    new = delta + v
                    return new
                elif once==True:
                    J_once = J(v0)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

            elif finite == True:
                if once ==False:
                    delta = np.linalg.solve(finite_difference(v, h), -F(v))
                    new = delta + v
                    return new
                elif once == True:
                    J_once = finite_difference(v0, h)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

        while k < max_iter and np.linalg.norm(xvals[-1]) < x_max and np.linalg.norm(F(xvals[-1])) > margin:
            try:
                new = solving_motor(xvals[-1])
                xvals.append(new)

            except sl.LinAlgError:
                y = xvals[-1] + step
                xvals.append(y)
            k += 1

        if np.linalg.norm(F(xvals[-1])) > margin or np.linalg.norm(xvals[-1]) > x_max:
            return (np.array((np.inf, np.inf)), k)

        else:
            return xvals[-1], k
    """
    def index_zero_motor(v,k, F=None, J=None, Array=[], Ind=[]):
        v_new   = Newt_omega(F, J, v)
        v_array = v_new[0]
        v_ind   = v_new[1]

        if len(Array) == 0 and len(Ind) == 0:
            if v_ind == max_iter:
                if (np.abs(F(v_array)) < margin).all():
                    Array.append(v_array)
                    Ind.append(v_ind)
                else:
                    pass

            elif v_ind != max_iter:
                Array.append(v_array)
                Ind.append(v_ind)
        else:
            for i, s in enumerate(Array):
                boolean = (np.abs(v_array -s) < margin).all()
                if boolean == True:
                    return v_ind
                else:
                    Array.append(v_array)
                    Ind.append(v_ind)
    """
    if plotting_method == "fml":
        def zero_collect_motor(F, J, v):
            v_a, v_i = Newt_motor(F, J, v)
            if (v_a==np.inf).any():
                return (np.inf, v_i)
            if not len(zero):
                zero.append(v_a)
                return (0,v_i)
            for ind, s in enumerate(zero):
                if np.linalg.norm(s-v_a) < margin:
                    return (ind, v_i)
            zero.append(v_a)
            return (len(zero)-1, v_i)


        A = np.zeros((N,N))
        L = tensor_motor(a,b,c,d,N)
        A = [[zero_collect_motor(F, J, L[i,j])[0] for i in range(N)] for j in range(N)]



    elif  plotting_method == "mine":
        A = np.zeros((N,N))
        L = tensor_motor(a,b,c,d,N)
        for i in range(N):
            for j in range(N):
                new = Newt_motor(F, J, L[j,i])[0]
                A[i,j] = np.abs(new[0]) + np.abs(new[1]**1)




    if PLOT == True:
        x = np.linspace(a,b,N)
        y = np.linspace(c,d,N)
        X,Y = np.meshgrid(x,y)
        fig, ax = plt.subplots()
        ax.pcolor(X, Y, A, cmap='plasma')

    elif PLOT == False:
        return A

    elif PLOT =="both":
        x = np.linspace(a,b,N)
        y = np.linspace(c,d,N)
        X,Y = np.meshgrid(x,y)
        fig, ax = plt.subplots()
        ax.pcolor(X, Y, A, cmap=f"{coloroption}")

        return zero


def Test_F(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([3*x**3-3*x*y**2 -1, 3*x**2*y -y**3])



def Test_J(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[9*x**2 - 3*y**2, -6*x*y],
                     [6*x*y, 3*x**2 - 3*y**2]])


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])



def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])

def Test_F2(X):
    X = np.array(X)
    return np.array([X[0]**3 - 3*X[0]*X[1]**2 - 2*X[0] -2, 3*X[0]**2*X[1] - X[1]**2 -2*X[1]])


def Test_J2(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2 -2 , -6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 2*X[1] -2]])


def Test_F3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([ x**8 -28*x**6*y**2 + 70*x**4*y**4 +15*x**4 -28*x**2*y**6 -90*x**2*y**6 + y**8 + 15*y**4 -16,  8*x**7*y -56*x**5*y**3 + 56*x**3*y - 8*x*y**7 -60*x*y**3])


def Test_J3(X):
    X = np.array(X)
    x = X[0]
    y = X[1]
    return np.array([[8*x**7 - 168*x**5*y**2 + 280*x**3*y**4 + 60*x**3 - 236*x*y**6, -56*x**6*y + 280*x**4*y**3 - 708*x**2*y**5 + 8*y**7 + 60*y**3 ],
                     [56*x**6*y - 280*x**4*y**3 + 168*x**2*y - 8*y**7 - 60*y**3, 8*x**7 - 168*x**5*y**2 + 56*x**3 - 56*x*y**6 - 180*x*y**2]])



#print(motor_integration(-2,2,-2, 2, 10, Test_F2, Test_J2, PLOT="both",plotting_method="mine", coloroption="Spectral"))
#print(motor_integration(-2,2,-2, 2,10, Test_F2, Test_J2, PLOT="both",plotting_method="fml", coloroption="gray"))
print(motor_integration(-2,2,-2,2,1000, Test_F, Test_J, PLOT="both",plotting_method="mine", coloroption="plasma"))
print(motor_integration(-2,2,-2, 2,1000, Test_F, Test_J, PLOT="both",plotting_method="fml", coloroption="plasma"))
#print(motor_integration(-2,2,-2, 2, 10, Test_F3, Test_J3, PLOT="both",plotting_method="mine", coloroption="inferno_r"))
#print(motor_integration(-2,2,-2, 2,10, Test_F3, Test_J3, PLOT="both",plotting_method="fml", coloroption="gray"))






















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










A = [1]
print(not  len(A))

def  collect_insert_motor(F, J, v,margin=1.e-4, max_iter=100, Array=[]):
    def Newt_motor(F, J, v0, max_iter=max_iter, margin=margin, step=step, h=1.e-9, finite=False, once=False):
        v0 = np.array(v0)
        xvals = [v0]
        V = xvals[-1]
        V_nrm = np.linalg.norm(xvals[-1])
        F_nrm = np.linalg.norm(F(xvals[-1]))
        x_max = 1/margin
        F_max = 1/margin
        k = 0

        def solving_motor(v):
            def finite_difference(v, h):
                J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                              [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                return J
            if finite == False:
                if once==False:
                    delta = np.linalg.solve(J(v), -F(v))
                    new = delta + v
                    return new
                elif once==True:
                    J_once = J(v0)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

            elif finite == True:
                if once ==False:
                    delta = np.linalg.solve(finite_difference(v, h), -F(v))
                    new = delta + v
                    return new
                elif once == True:
                    J_once = finite_difference(v0, h)
                    delta = np.linalg.solve(J_once, -F(v))
                    new = delta + v
                    return new

        while k < max_iter and np.linalg.norm(xvals[-1]) < x_max and np.linalg.norm(F(xvals[-1])) > margin:
            try:
                new = solving_motor(xvals[-1])
                xvals.append(new)

            except sl.LinAlgError:
                y = xvals[-1] + step
                xvals.append(y)
            k += 1

        if np.linalg.norm(F(xvals[-1])) > margin or np.linalg.norm(xvals[-1]) > x_max:
            return (np.array((np.inf, np.inf)), k)

        else:
            return xvals[-1], k

    x_max = 1/margin
    if v.ndim == 1:
        new = Newt_motor(F,J,v)
        v_array = new[0]
        v_ind   = new[1]
        if np.linalg.norm(F(v_array)) < margin:
            return v_array
        elif np.linalg.norm(F(v_array)) > margin or np.linalg.norm(v_array) > x_max:
            raise Exception(f"{v} does not converge to any zero")

    elif v.ndim == 3:

        def for_motor(A):
            matrix = np.array((A.shape[0], A.shape[1]))
            for i in range(A.shape[0]):
                for j in range(A.shape[1]):
                    new = Newt_motor(F, J, A[i,j])
                    n_abs = np.linalg.norm(new[0])
                    for s in Array:
                        if np.abs(s - n_abs) < margin:
                            matrix[i,j] = x_max
                        else:
                            matrix[i,j] = n_abs
            return matrix

        def gauge_norm_motor(vector):
                first = Newt_motor(F, J, vector)
                f_a = first[0]
                f_i = first[1]
                if np.linalg.norm(F(f_a)) < margin:
                    return True, f_a
                else:
                    return False
        var1 = gauge_norm_motor(v[0,0])
        if var1[0] == True:
            f_abs = np.linalg.norm(var[1])
            Array.append(f_abs)
        else:
            for i in range(1, v.shape[0]):
                for j in range(1, v.shape[1]):
                    var2 = gauge_norm_motor(v[i,j])
                    if var2[0] ==True:
                        f_abs = np.linalg.norm(var[1])
                        Array.append(f_abs)

            if len(Array) == 0:
                raise Exception("the meshgrid contains no zero at all")

        M = for_motor(v)
        return M


def Test_F1(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])


def Test_J1(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])


print(collect_insert_motor(Test_F1, Test_J1, np.array([3., 4])))























def index_zero_motor(v,k, F=None, J=None, Array=[], Ind=[]):
    v_new   = Newt_omega(F, J, v)
    v_array = v_new[0]
    v_ind   = v_new[1]

    if len(Array) == 0 and len(Ind) == 0:
        if v_ind == max_iter:
            if (np.abs(F(v_array)) < margin).all():
                Array.append(v_array)
                Ind.append(v_ind)
            else:
                pass

        elif v_ind != max_iter:
            Array.append(v_array)
            Ind.append(v_ind)
    else:
        for i, s in enumerate(Array):
            boolean = (np.abs(v_array -s) < margin).all()
            if boolean == True:
                return v_ind
            else:
                Array.append(v_array)
                Ind.append(v_ind)



"""
Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap,
 CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r,
 Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r,
 PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r,
 RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1,
 Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r,
 YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone,
 bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r,
 cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat,
 gist_heat_r,
 gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg,
 gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r,
 inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean,
 ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic,
 seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b,
 tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
"""

u

"""
Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap,
 CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r,
 Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r,
 PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r,
 RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1,
 Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r,
 YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone,
 bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r,
 cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat,
 gist_heat_r,
 gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg,
 gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r,
 inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean,
 ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic,
 seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b,
 tab20b_r, tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r, winter, winter_r
"""

x, y = symbols('x y', real=True)
f = y*(x/(1+x)) +1 +x

d_x = diff(f, x)
d_y = diff(f, y)

print(d_x)
print(d_y)


def g(x):
    return x**2

f = lambdify((x,y),g(x), "numpy")
print(f(1, 2))






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











def LOL(F, J, v0, max_iter=10, margin=1.e-4, step=0.001, h=1.e-9, approx=False, once=False):
    v0 = np.array(v0)
    xvals = [v0]
    V = xvals[-1]
    V_nrm = np.linalg.norm(xvals[-1])
    F_nrm = np.linalg.norm(F(xvals[-1]))
    x_max = 1/margin
    F_max = 1/margin
    k = 0

    def Motor_1(v):
        def finite_difference(v, h):
            J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                      [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
            return J
        if approx == False:
            if once==False:
                delta = np.linalg.solve(J(v), -F(v))
                new = delta + v
                return new
            elif once==True:
                J_once = J(v0)
                delta = np.linalg.solve(J_once, -F(v))
                new = delta + v
                return new

        elif approx == True:
            if once ==False:
                delta = np.linalg.solve(finite_difference(v, h), -F(v))
                new = delta + v
                return new
            elif once == True:
                J_once = finite_difference(v0, h)
                delta = np.linalg.solve(J_once, -F(v))
                new = delta + v
                return new



    while k < max_iter and np.linalg.norm(xvals[-1]) < x_max or np.linalg.norm(F(xvals[-1])) < F_max:
        try:
            new = Motor_1(xvals[-1])
            xvals.append(new)

        except sl.LinAlgError:
            y = xvals[-1] + step
            xvals.append(y)


        k += 1

 #   if np.linalg.norm(F(xvals[-1])) > F_max or np.linalg.norm(xvals[-1]) > x_max:
  #      return (np.array((np.inf, np.inf)), k)

    else:
        return xvals[-1], k



def Test_F(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])

def Test_J(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])


c = np.array([0.000001,0.0000001])
print(LOL(Test_F,Test_J, c))







def LEL(F, J, v0):
    v0 = np.array(v0)
    xvals = [v0]
    V = xvals[-1]
    V_nrm = np.linalg.norm(xvals[-1])
    F_nrm = np.linalg.norm(F(V))
    x_max = 100
    k = 0

    while k < 10 or np.linalg.norm(xvals[-1]) < x_max or np.linalg.norm(F(xvals[-1])) > 0.0001 :
        try:
            delta = np.linalg.solve(J(xvals[-1]), -F(xvals[-1]))
            new = delta + xvals[-1]
            xvals.append(new)
            k += 1
        except sl.LinAlgError:
            y = xvals[-1] + 0.001
            xvals.append(y)

    return xvals


def Test_F(X):
    X = np.array(X)
    return np.array([X[0]**3 -3*X[0]*X[1]**2 -1, 3*X[0]*X[1] - X[1]**3])


def Test_J(X):
    X = np.array(X)
    return np.array([[3*X[0]**2 - 3*X[1]**2, 6*X[0]*X[1]],
                     [6*X[0]*X[1], 3*X[0]**2 - 3*X[1]**2]])



c = np.array([1.,8])
print(LEL(Test_F,Test_J, c))





def N_function_omega(F, x0, J=None, max_iter=10000000, margin=1.e-4, step=0.001, h=1.e-9, standard=False):
    xvals =[x0]
    X_L = xvals[-1]
    X_nrm = np.linalg.norm(X_L)
    F_nrm = np.linalg.norm(F(*X_L))
    zeros = []
    x_max = 1/margin
    #J_once = J(*x0)


    def finite_difference(x0, h):
        J = np.array([[(F[0](x0[0]+h,x0[1]) -F[0](*x0))/h ,    (F[0](x0[0],x0[1]+h) -F[0](*x0))/h],
                      [(F[1](x0[0]+h,x0[1]) -F[1](*x0))/h ,    (F[1](x0[0],x0[1]+h) -F[1](*x0))/h]])
        return J
    if not J:
            J = lambda x0: finite_difference(x0, h)

    def V_to_new(V):
        if standard==True:
            delta = np.linalg.solve(J(*V), -F(*V))
            x_new = delta + V
            return x_new

        elif standard==False:
            delta = np.linalg.solve(J_once, -F(*V))
            x_new = delta + V
            return x_new

        elif not J:
            J = lambda x0: finite_difference(x0, h)
            delta = finite_difference(V, h)
            x_new = delta + V
            return x_new


    k = 0
    while k < max_iter and X_nrm < x_max and F_nrm < margin:
        try:
            new = V_to_new(X_L)
            xvals.append(new)

        except sl.LinAlgError:
            y = X_L + step
            xvals.append(X_L)

         #if this line's reached <=> J(*X_L) was not singular, x_L is immediately stored
        k += 1

    if F_nrm > margin or X_nrm > x_max:
        return (np.array((np.inf, np.inf)), max_iter)
    return (X_L, k)

















def finite_difference(self, x0, h):
    F = self.F(*x0)
    J = np.array([[(F[0](x0[0]+h,x0[1]) -F[0](*x0))/h ,    (F[0](x0[0],x0[1]+h) -F[0](*x0))/h],
                  [(F[1](x0[0]+h,x0[1]) -F[1](*x0))/h ,    (F[1](x0[0],x0[1]+h) -F[1](*x0))/h]])
    return J




def Newton_omega(self,x0, max_iter=100, margin=1.e-4, step=0.001, standard=True):
    x0 = np.array(x0)
    xvals =[x0]
    X_L = xvals[-1]
    X_nrm = np.linalg.norm(X_L)
    F_nrm = np.linalg.norm(self.F(*X_L))
    zeros = []
    x_max = 1/margin
    J_once = self.J(*x0)


    def V_to_new(V):
        if standard==True:
            delta = np.linalg.solve(self.J(*V), -self.F(*V))
            x_new = delta + V
            return x_new

        elif standard==False:
            delta = np.linalg.solve(J_once, -self.F(*V))
            x_new = delta + V
    k = 0
    while k < max_iter and X_nrm < x_max and F_nrm < margin:
        try:
            new = V_to_new(X_L)
        except sl.LinAlgError:
            X_L = X_L + step

        xvals.append(X_L) #if this line's reached <=> self.J(*X_L) was not singular, x_L is immediately stored
        k += 1

    if f_nrm > margin or X_nrm > x_max:
        return (np.array((np.inf, np.inf)), max_iter)

    return (X_L, k)









A = np.array([[[[0.,  0. ]
               [0.,  0.5]
               [0.,  1. ]],

               [[0.5, 0. ]
                [0.5, 0.5]
                [0.5, 1. ]],

                [[1.,  0. ]
                 [1.,  0.5]
                 [1.,  1. ]]]])









def plot_array(a,b,c,d,N):

    def make_tensor(a,b,c,d,N):
        """
        This function creates a block Tensor, e.g. a tensor of shape (i,j,k).
        Thsi is to make a matrix P tha has tuples as matrix elemts, that will be set into
        the Newtons method.

        PARAMETERS:
            a,b,c,d,N (float or int):
                are the interval parameters that create a rectangle of dimesion [a,b]x[c,d]
                that defines the domain that the newtin method will operate over.

        RETURNS:
            L (array):
                the shape of this array is (N, N, 2)
                and the array L[i,j] = p[i,j] = (Y[i,j], X[i,j])  and can be set into the newtin method.
        """
        x = np.linspace(a,b,N)
        y = np.linspace(c,d,N)
        L = []
        X,Y = np.meshgrid(x,y)

        for i in range(N):
            L.append([])

        for i in range(X.shape[0]):
            for j in range(X.shape[1]):
                L[i].append((Y[i,j], X[i,j]))
        L = np.array(L)
        return  L


    def N_function_omega(F, x0, J=None, max_iter=10000000, margin=1.e-4, step=0.001, h=1.e-9, standard=False):
        xvals =[x0]
        X_L = xvals[-1]
        X_nrm = np.linalg.norm(X_L)
        F_nrm = np.linalg.norm(F(*X_L))
        zeros = []
        x_max = 1/margin
        #J_once = J(*x0)


        def finite_difference(x0, h):
            J = np.array([[(F[0](x0[0]+h,x0[1]) -F[0](*x0))/h ,    (F[0](x0[0],x0[1]+h) -F[0](*x0))/h],
                           [(F[1](x0[0]+h,x0[1]) -F[1](*x0))/h ,    (F[1](x0[0],x0[1]+h) -F[1](*x0))/h]])
            return J
        if not J:
                J = lambda x0: finite_difference(x0, h)

        def V_to_new(V):
            if standard==True:
                delta = np.linalg.solve(J(*V), -F(*V))
                x_new = delta + V
                return x_new

            elif standard==False:
                delta = np.linalg.solve(J_once, -F(*V))
                x_new = delta + V
                return x_new

            elif not J:
                J = lambda x0: finite_difference(x0, h)
                delta = finite_difference(V, h)
                x_new = delta + V
                return x_new
        k = 0
        while k < max_iter and X_nrm < x_max and F_nrm < margin:
            try:
                new = V_to_new(X_L)
            except sl.LinAlgError:
                X_L = X_L + step

            xvals.append(X_L) #if this line's reached <=> J(*X_L) was not singular, x_L is immediately stored
            k += 1

        if F_nrm > margin or X_nrm > x_max:
            return (np.array((np.inf, np.inf)), max_iter)
        return (X_L, k)

    def zero_index(V, zero=[], x_iter_list=[]):
        new = Newton(F, J, V)
        l = len(zero)

        if l==0 and new != np.array((np.inf, np.inf)):
            zero.append(new[0])
            x_iter_list.append(new[1])
        elif new == np.array((np.inf, np.inf)):
            return np.inf
        else:
            for ind, s in enumerate(zero):
                if (np.abs(s - new) < margin).all() == True:
                    return ind
                else:
                    zero.append(new[0])
                    x_iter_list.append(new[1])

    A = np.zeros((N,N))
    L = make_tensor(a,b,c,d,N)
    for i in range(N):
        for j in range(N):
            new = Newton2(F, J, L[i,j])
            if isinstance(new, float):
                A[i,j] = new

    fig, ax = plt.subplots()
    ax.pcolor(X, Y, A, cmap='plasma')
























L= make_array(0,1,0,1,3)

print(L[0,1])
print(L)
print("------------------")







def parallell(x,y):
    X,Y = np.meshgrid(x,y)
    list_ = np.zeros((X.shape))
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            a = np.array([Y[i,j], X[i,j]])
            list_[i,j] = a
    return list_

print(parallell(Q,E))











def Taylor_Pol(x,a,n, func=None, alpha=None):
    T = x-a

    def fac(s):
        if s == 0:
            return 1
        else:
            return s*fac(s-1)

    if func == "sin":
        M = [T**(2*k +1)/fac(2*k +1)*(-1)**(k) for k in range(n)]
        return sum(M)

    if func == "cos":
        M = [T**(2*k)/fac(2*k)*(-1)**(k) for k in range(n)]
        return sum(M)

    if func == "cos":
        M = [T**(k)/fac(k) for k in range(n)]
        return sum(M)

    if func == "arctan":
        M = [T**(2*k +1)/(2*k +1)*(-1)**k for k in range(n)]
        return sum(M)

    if func == "log":
        M = [T**(k +1)/(k +1)*(-1)**k for k in range(n)]
        return sum(M)

    if func == "polynomial":
        M = [fac(alpha)/fac(alpha -k)*T**(k)/fac(k) for k in range(n)]
        return sum(M)






print(Taylor_Pol(np.pi/3,4, func ="cos"))

























def Deriv(y, n):
    a = sp.symbols('a', real=True)
    list_ = [cos(a)]
    f = list_[-1]
    for i in range(n):
        d_a = diff(f, a)
        list_.append(d_a)
    fp = sp.lambdify(a, list_[-1])
    a = y
    return fp(a)

print(Deriv(pi/3, 3))

def Taylor_expansion(x, deg, interpolation_point):

    def Deriv(y, n):
        a = sp.symbols('a', real=True)
        list_ = [cos(a)]
        f = list_[-1]
        for i in range(n):
            d_a = diff(f, a)
            list_.append(d_a)
        fp = sp.lambdify(a, list_[-1])
        a = y
        return fp(a)

    def factorial(n):
        z =1
        if n ==0:
            return 1
        else:
            for i in range(n):
                z *= i
            return z
    List = []
    for k in range(deg):
        val = Deriv(interpolation_point,k)*((x-interpolation_point)**k)/factorial(k)
        List.append(val)
    return List
print(Taylor_expansion(0.01, 4, 0))






def fp(x):
    partial_deriv(1)
fp = lambdify((x), fp(x))

print(fp(1))
print(partial_deriv(1))




def ODE_solver(f, fp, t0, u0, te, N=None, h = 4):
    U = [[u0], [t0]]
    Y = [[u0], [t0]]
    delta_t = np.abs(te-t0)

    while np.abs(U[1][0] - U[1][-1]) <= delta_t:
        RadiusOfCurvature = fp(U[0][-1], U[1][-1])
        u_new = U[0][-1] + h*RadiusOfCurvature*f(U[0][-1], U[1][-1])
        t_new = U[1][-1] + h*RadiusOfCurvature

        U[0].append(u_new)
        U[1].append(t_new)

    while np.abs(Y[1][0] - Y[1][-1]) <= delta_t:
        u_new = Y[0][-1] + h*f(Y[0][-1], Y[1][-1])
        t_new = Y[1][-1] + h

        Y[0].append(u_new)
        Y[1].append(t_new)

    def y(x):
        return -3 -x - 2/(1+x) +5*np.exp(x)/(1+x)
    plt.figure()
    plt.plot(U[1], U[0], label = "curvature regulated euler",color =  "r")
    plt.plot(Y[1], Y[0], label = "standard euler",  color = "b")
    x = np.linspace(t0, te, (len(Y[1]) + len(U[1]))/2)
    plt.plot(x, y(x), label = "actual solution", color = "g")
    plt.title("comparing solutions to ODE")
    plt.legend(loc = "best")


def f(y,x):
    return y*(x/(1+x)) +1 +x

def fp(y,x):
    return  np.abs(x/(x + 1)*(y*(x/(1+x)) +1 +x) -x*y/(x + 1)**2 + y/(x + 1) + 1)/(1+ (y*(x/(1+x)) +1 +x)**2)

ODE_solver(f,fp, 0, 0, 10)


class TaylorPol:
    def __init__(self, f, a, *coeff):
        self.coeff = coeff
        self. f = f
        self.a = a

    def __calc_coeff__(self):
        f = self.f












def PLOOT(l,u, N, delta, p=None):

    if not p:
        p = lambda I: 3*I**3 -2*I**2 - 5*I -1

    def PLOT(l, u, N, delta):
        vp = np.vectorize(p)
        xl = np.linspace(l,u,N)
        xu = np.linspace(l,u,N) + delta
        Tl = vp(xl)
        Tu = vp(xu)

        plt.figure()
        plt.plot(xl, Tl, label = "transformed lower bounds", color ="g")
        plt.plot(xl, Tu, label = "transformed upper bounds", color ="r")
        plt.xlabel("x")
        plt.ylabel("p(I)")
        plt.title("$p(I) =  3I^3 -2I^2 - 5I -1$, I=interval(x,x+{delta})")
        plt.xlim(0,1)
        plt.ylim(-10,4)

    if not l and not u and not N and not delta:
        l = 0
        u = 1
        N = 1000
        delta = 0.5
        return PLOT(l,u,N,delta)
    else:
        return PLOT(l,u,N,delta)

PLOOT(l=None, u=None, N=None, delta=None)





M = np.array([[1,2.],
              [1,2]])
b = np.array([4,5])
try:
    c =  np.linalg.solve(M, b)
except sl.LinAlgError:
     print('it works')

fig , ax = plt.subplots(figsize = (6, 3))
x = np.linspace(-2*np.pi, 2*np.pi, 100)
ax.plot(x, np.sin(x), "r--", label = "sine")
ax.plot(x, np.cos(x), "b", label = "cosine")
print("The  number  of line  plots  is", len(ax.lines ))
print("Each  line  plot is of the  type", type(ax.lines[0]))

"""
LECTURE 9?
HOW TO PLOT IN AN OBJECT ORIENTED WAY:
    In this unit we study the different components of matplotlib and see how we can plot in
    an object oriented way.

    THrer are two steps:
        - contruction of the "internal" representation of a plt.
        - representing a lplot  on the screen in  a window or the filesystem as a file in a special graphical
          fileforamt.
    THE TWO PROGRAMMING PARADIGMS:
        function driven & object driven.
        the subplot commmand example: subplot(ijk) the ij gives the simension of the "maytrox" of subplots.
        aka, subplot(32k) would then mean that we have a matrix with 3 rows and 2 columns. therefore (not zero-based) 1<= k <= 2.

        to create one figure we can just write:
            fig, ax = plt.subplot()
            # ax is the literal name of the window,
        to create mutiple  small windows we instead write:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows = 2, ncols = 2)
            ax1 , ax2 ax3 and ax4 ae the names of teh respective subplots

        THE FOLLOWING PLOTS ARE THE SAME:

        USING FUNCTIONS:
        x = linspace(-2*pi, 2*pi, 100)figure ()# not  needed  when  using  Spyder
        plt.plot(x, sin(x/2), label = ’sine’)
        plt.plot(x, cos(x/2), label = ’$\cos(x)$’)
        plt.legend(loc = ’lower  center ’, fontsize = ’small’)

        USING METHODS
        x = linspace(-2*pi, 2*pi, 100)
        fig , ax = plt.subplots()
        ax.plot(x, sin(x/2), label = ’sine’)
        ax.plot(x, cos(x/2), label = ’$\cos(x)$’)
        ax.legend(loc = ’lower  center ’, fontsize = ’small’)

        matplotlib has a plot object typ with its own methods. THis is therefore a class.
        one method is a "save" method where you can save a plot as a file.





    method plotting first we import matplotlib:
        "import matplotlib.pyplot as plt"


"""
def Newton3(F, J, x0, max_iter=100, margin=1.e-4, step=0.001):
    x0 = np.array(x0)
    xvals =[x0]
    X_L = xvals[-1]
    X_nrm = np.linalg.norm(X_L)
    F_nrm = np.linalg.norm(F(*X_L))
    zeros = []
    x_max = 1/margin

    def solve(V):
        delta = np.linalg.solve(J(*V),-F(*V))
        X_n = delta + V
        xvals.append(X_n)














def Newton2(F, J, x0, max_iter=100, margin=1.e-4, step=0.001):
        x0 = np.array(x0)
        xvals =[x0]
        X_L = xvals[-1]
        X_nrm = np.linalg.norm(X_L)
        F_nrm = np.linalg.norm(F(*X_L))
        zeros = []
        x_max = 1/margin

        def V_to_new(V):
            new = np.linalg.solve(J(*V), -F(*V))
            return new
        k = 0
        while k < max_iter or X_nrm < x_max or F_nrm < margin:
            try:
                delta = V_to_new(X_L)
                xvals.append(new)
            except sl.LinAlgError:
                X_L = X_L + step
            k +=1
        if np.abs(F_nrm) > tol or np.abs(X_nrm) > x_max:
            return (np.array((np.inf, np.inf)), max_iter)

        return (X_L, k)

def zero_index(x0, zero=[], x_iter_list=[]):
    new = Newton(F, J, x0)
    l = len(zero)

    if l==0 and new != np.array((np.inf, np.inf)):
        zero.append(new[0])
        x_iter_list.append(new[1])
    elif new == np.array((np.inf, np.inf)):
        return np.inf
    else:
        for s in zero:
            if (np.abs(s - new) < margin).all() == True:
                ind = zero.index(s)
                return ind
            else:
                zero.append(new[0])
                x_iter_list.append(new[1])



def Newton3(F, J, x0, max_iter=100, margin=1.e-4, step=0.001):
        x0 = np.array(x0)
        xvals =[x0]
        X_L = xvals[-1]
        X_nrm = np.linalg.norm(X_L)
        F_nrm = np.linalg.norm(F(*X_L))
        zeros = []
        x_max = 1/margin

        def V_to_n(V):
            N = np.matmul(J(*V),V) -F(*V)
            n = np.linalg.solve(J, N)
            return n
        k = 0
        while k < max_iter or X_nrm < x_max or F_nrm < margin:
            try:
                new = V_to_n(X_L)
                xvals.append(new)
            except LinAlgError:
                np.abs(np.linalg.det(J(*X_L)))
                y = X_L + step
                if np.abs(np.linalg.det(J(*y))) != 0.:
                    new = V_to_n(y)
                    xvasl.append(new)
                else:
                    while not np.abs(np.linalg.det(J(*y))):
                        y += step
                    new = V_to_N(y)
                    xvals.append(y)

            k +=1
        if np.abs(F_nrm) > tol or np.abs(X_nrm) > x_max:
            return np.array((np.inf, np.inf))

        return (X_L, k)

def zero_index(x0, zero=[], x_iter_list=[]):
    new = Newton(F, J, x0)
    l = len(zero)

    if l==0 and new != np.array((np.inf, np.inf)):
        zero.append(new[0])
        x_iter_list.append(new[1])
    elif new == np.array((np.inf, np.inf)):
        return np.inf
    else:
        for ind, s in enumerate(zero):
            if (s - new < margin).all() == True:
                ind = zero.index(s)
                return ind
            else:
                zero.append(new[0])
                x_iter_list.append(new[1])


























class Fractal2D:
    def __init__(self, F, J=None, h=0.001):
        """
        Task 1
        """
        self.F = F
        self.J = J
        if not J:
            self.J = lambda x,y: self.finite_difference(x,y,h)

        self.zeros =[]
        self.x_iter_list =[]

    def Newton(self, x0, max_iter = 1000, margin=1.e-4, step =0.001):
        """
        Task 2
        """
        x0 = np.array(x0)
        xvals =[x0]
        X_L = xvals[-1]
        X_nrm = np.linalg.norm(X_L)
        F_nrm = np.linalg.norm(self.F(*X_L))
        zeros = []
        x_max = 1/margin

        def V_to_n(V):
            N = np.matmul(self.J(*V) -self.F(*V))
            return N
        k = 0
        while k < max_iter or X_nrm < x_max or F_nrm < margin:
            try:
                new = V_to_n(X_L)
                xvals.append(new)
            except LinAlgError:
                y = X_L + step
            k +=1

        if np.abs(F_nrm) > tol or np.abs(X_nrm) > x_max:
                return np.array((np.inf, np.inf), max_iter)

        return (X_L, k)



    def zero_index(self, x0, margin=1.e-4):
        """
        Task 3
        """
        new = Newton(self, x0)
        l = len(self.zero)

        if l==0 and new[0] != np.array((np.inf, np.inf)):
            self.zero.append(new[0])
            self.x_iter_list.append(new[1])
        elif new == np.array((np.inf, np.inf)):
            return np.inf
        else:
            for ind, s in enumerate(self.zero):
                if (np.abs(s - new[0]) < margin).all() == True:
                    #ind = self.zero.index(s)
                    return ind
                else:
                    self.zero.append(new[0])
                    self.x_iter_list.append(new[1])



    def plot_array(self,a,b,c,d,N):

        def make_tensor(a,b,c,d,N):
            """
            This function creates a block Tensor, e.g. a tensor of shape (i,j,k).
            Thsi is to make a matrix P tha has tuples as matrix elemts, that will be set into
            the Newtons method.

            PARAMETERS:
                a,b,c,d,N (float or int):
                    are the interval parameters that create a rectangle of dimesion [a,b]x[c,d]
                    that defines the domain that the newtin method will operate over.

            RETURNS:
                L (array):
                    the shape of this array is (N, N, 2)
                    and the array L[i,j] = p[i,j] = (Y[i,j], X[i,j])  and can be set into the newtin method.
            """
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            L = []
            X,Y = np.meshgrid(x,y)
            for i in range(N):
                L.append([])

            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    L[i].append((Y[i,j], X[i,j]))
            L = np.array(L)
            return  L


        def Newton2(self,x0, max_iter=100, margin=1.e-4, step=0.001):
            x0 = np.array(x0)
            xvals =[x0]
            X_L = xvals[-1]
            X_nrm = np.linalg.norm(X_L)
            F_nrm = np.linalg.norm(self.F(*X_L))
            zeros = []
            x_max = 1/margin

            def V_to_new(V):
                delta = np.linalg.solve(self.J(*V), -self.F(*V))
                x_new = delta + V
                return x_new
            k = 0
            while k < max_iter and X_nrm < x_max and F_nrm < margin:
                try:
                    new = V_to_new(X_L)
                    xvals.append(new)
                except sl.LinAlgError:
                    X_L = X_L + step
                k +=1
            if np.abs(F_nrm) > tol or np.abs(X_nrm) > x_max:
                return (np.array((np.inf, np.inf)), max_iter)

            return (X_L, k)

        def zero_index(V, zero=[], x_iter_list=[]):
            new = Newton(F, J, V)
            l = len(zero)

            if l==0 and new != np.array((np.inf, np.inf)):
                zero.append(new[0])
                x_iter_list.append(new[1])
            elif new == np.array((np.inf, np.inf)):
                return np.inf
            else:
                for ind, s in enumerate(zero):
                    if (np.abs(s - new) < margin).all() == True:
                        return ind
                    else:
                        zero.append(new[0])
                        x_iter_list.append(new[1])

        A = np.zeros((N,N))
        L = make_tensor(a,b,c,d,N)
        for i in range(N):
            for j in range(N):
                new = Newton2(F, J, L[i,j])
                if isinstance(new, float):
                    A[i,j] = new

        fig, ax = plt.subplots()
        ax.pcolor(X, Y, A, cmap='plasma')





















def newton(f,fp, x0, max_iter =400, div_deterrent =10**(10), margin =1.e-9,slope_near_zero = 1.e-9, step = 10**9, climb_steps= 10**(4)):
    x =[x0]
    def check_slope(fp, X, Where_):
        SLOPE = fp(X)
        if np.abs(SLOPE) < slope_near_zero:
            for j in range(climb_steps):
                X += step
                Yp = fp(X)
                if Yp > slope_near_zero:
                    x[-1] = X
                    return x


            if Where_ =="start of code":
                raise Exception("slope near zero at initial evaluation")

            elif Where_ =="middle of code":
                raise Exception("slope near zero at later evaluations")

    check_slope(fp, x0, Where_="start_of_code")


    for i in range(max_iter):
        x_new = x[-1] -f(x[-1])/fp(x[-1])
        x.append(x_new)
        delta = np.abs(x[-1] - x[-2])
        travel_len = np.abs(x[0] - x[-1])
        check_slope(fp, x[-1], Where_="middle of code")


        if travel_len < div_deterrent:
            if i < max_iter:
                if delta < margin:
                    return x[-1], f"convergent, {x[-1]} = x[-1] < margin"
            elif i == max_iter:
                if delta > margin:
                    return  f"divergent, {x[-1]} = x[-1] > margin"
        elif travel_len >= div_deterrent:
            raise Exception(f"x-values are travelling too far away {x[-1]}")

    return x[-1]






class Nameless:
    def __init__(self, b, W, margin=1.e-6):
        """
        b should be a matrix whoese rows are the basis vectors of
        a vectorspace.

        W should be a matrix whose rows contain the apporeriate scalars
        such that when one of the rows is multipied with the matrix b,
        one should then get a matrix whose rwos, when summed, give you a
        vector that is a linear combination of the basis vectors.
        """

        if not isinstance(b, np.ndarray):
            raise TypeError("b should be a matrix")
        elif np.abs(ls.det(d)) < margin:
            raise Exception("basis vetcors are almost linearly dependent")

        if not isinstance(W, np.ndarray):
            raise TypeError("W should be a matrix")


        VectorMatrices =[]
        L = W.shape[1]
        for i in range(1, L):
            V = b *W[i]
            VectorMatrices.append(V)

        def Functional(b):
            if b != np.eye(L):
                raise Exception("b is not the identity matrix")
            else:
                Delta = b
                return Delta

        self.KroneckerDelta = Functional(b)

    def __CoordTransform__(self, X, XP, h=1.e-10):
        """
        This method is supposed to take the identity matrix of in tensor analysis
        mostly refered to as the "Kronecker delta" and transfrom it


        """

        def phi_psi(X,Xp):
            list_ =[]
            for x in X:
                x = x(*Xp)
                list_append(x)
            Epsilon = np.array(list_)
            return list_

        def fin_diff(f, x):
            """
            this function computes the derivative at a point x using finite difference.
            """
            fp = (f(x+h) -f(x))/h
            return fp

        def many_dot(reference, *M):
            """
            This function matrix-multplies the matrices in the list M

            PARAMETERS:
                reference (array):
                    this matrix is used to make sure that when we compare the elements
                    in M with it, the function raises an exception if one of the
                    elements in M are not the same dimension, by using M[0] as the reference.

                M (list):
                    this list contains the matrices that you want to multiply.

            RETURNS:
                Var[-1] (array):
                    Var[-1] is the resultant matrix.

            NOTICE:
                the order of the elements in M matter, matrix multiplication is genreally not
                commutative and this function runs in accordance with these rules.
            """
            if type(reference) != np.ndarray:
                raise ValueError(f"the reference, {reference} should be a an array")
            elif np.ndim(reference) != 2:
                raise ValueError("this function only handles matrices")
            else:
                Var =[np.eye(reference.shape[0])]
                for s in M:
                    if s.shape != reference.shape:
                        raise Exception(f"{M} contains a matrix whose shape differs from the reference")
                    new =np.matmul(Var[-1],s)
                    Var.append(new)
                return Var[-1]

        phi = phi_psi(X, Xp)
        L = len(phi)
        Jacobian = np.zeros((L,L))

        for x in phi:
            empty =[]
            for xp in Xp:
                x_diff = fin_diff(x, xp, h)
                empty.append(x_diff)
            r = np.array(empty)
            ind = phi.index(x)
            J[ind] = r


        M = [J.T, self.KroneckerDelta, J]
        G = many_dot(M[0], *M)
        self.G = G
        return Nameless(self.G)

np.linalg.solve
np.linalg.LinAlgError































"""
x, y = symbols('x y', real=True)
f = y*(x/(1+x)) +1 +x
d_x = diff(f, x)
d_y = diff(f, y)

print(d_x)
print(d_y)
"""



def ODE_solver(f, fp, t0, u0, te, N=None, h = 4):
    U = [[u0], [t0]]
    Y = [[u0], [t0]]
    delta_t = np.abs(te-t0)

    while np.abs(U[1][0] - U[1][-1]) <= delta_t:
        RadiusOfCurvature = fp(U[0][-1], U[1][-1])
        u_new = U[0][-1] + h*RadiusOfCurvature*f(U[0][-1], U[1][-1])
        t_new = U[1][-1] + h*RadiusOfCurvature

        U[0].append(u_new)
        U[1].append(t_new)

    while np.abs(Y[1][0] - Y[1][-1]) <= delta_t:
        u_new = Y[0][-1] + h*f(Y[0][-1], Y[1][-1])
        t_new = Y[1][-1] + h

        Y[0].append(u_new)
        Y[1].append(t_new)

    def y(x):
        return -3 -x - 2/(1+x) +5*np.exp(x)/(1+x)
    plt.figure()
    plt.plot(U[1], U[0], label = "curvature regulated euler",color =  "r")
    plt.plot(Y[1], Y[0], label = "standard euler",  color = "b")
    x = np.linspace(t0, te, (len(Y[1]) + len(U[1]))/2)
    plt.plot(x, y(x), label = "actual solution", color = "g")
    plt.title("comparing solutions to ODE")
    plt.legend(loc = "best")


def f(y,x):
    return y*(x/(1+x)) +1 +x

def fp(y,x):
    return  np.abs(x/(x + 1)*(y*(x/(1+x)) +1 +x) -x*y/(x + 1)**2 + y/(x + 1) + 1)/(1+ (y*(x/(1+x)) +1 +x)**2)

ODE_solver(f,fp, 0, 0, 10)











class IVP(self, f, t0, u0, te, N, h, fp=None):
    self.f = f
    self.u0 = u0
    self.interval = [t0, te]
    if not fp:
        self.grid = np.linspace(t0, te, N)
        self.interval = (te-t0)/N
    else:
        pass











def many_dot(reference, *M):
    if type(reference) != np.ndarray:
        raise ValueError(f"the reference, {reference} should be a an array")
    elif np.ndim(reference) != 2:
        raise ValueError("this function only handles matrices")
    else:
        Var =[np.eye(reference.shape[0])]
        for s in M:
            if s.shape != reference.shape:
                raise Exception(f"{M} contains a matrix")
            new =np.matmul(Var[-1],s)
            Var.append(new)
        return Var[-1]

m = np.ones((3,3))
L =[m,m,m]
z = many_dot(L[0], *L)
print(z)









m = np.ones((3,3))
z = np.reduce(numpy.dot, [m,m,m])
print(z)


A =[1,2,3]
B =[1,1,1]

for x in A:
    empty=[]
    for y in B:
        c = x*y
        empty.append(c)
    r = np.array(empty)
    ind = A.index(x)
    m[ind] = r
print(m)










L =[x,y,z]
B =[]
for s in L:
    s = str(s)
    B.append(s)

Q = " ".join(B)
B = symbols(Q, real=True)
f = 4*x*y + x*sin(z) + x**3 + z**8*y
E = diff(f, B[0])
E = str(E)
E = Symbol(E)
print(type(E))


def f(*L):
    L =[x,y,z]
    B =[]
    for s in L:
        s = str(s)
        B.append(s)

    Q = " ".join(B)
    B = symbols(Q, real=True)
    f = 4*x*y + x*sin(z) + x**3 + z**8*y
    E = diff(f, B[0])
    E = str(E)
    E = Symbol(E)

    return E
print(f(1,1,1))







x, y, z = symbols('x y z', real=True)
f = 4*x*y + x*np.sin(z) + x**3 + z**8*y
E = diff(f, x)
print(E)

list1 = ['1', '2', '3']
str1 = ' '.join(list1)



L =["x","y","z"]
Q = " ".join(L)


L = symbols(Q, real=True)
f = 4*x*y + x*sin(z) + x**3 + z**8*y
E = diff(f, L[0])
print(E)


L =["x","y","z"]
*L = sp.symbol("*L")
print(sp.diff(x*y*z))



C = []
A = np.ones((3,3))
B = np.array([[1,2,3],
              [4,5,6],
              [7,8,9]])
print(B*B)

E = A*B[:1].reshape(-1,1)
print(E)
print(E[3])
Q = np.ones((3,8))
print(Q.shape)


e = "yess"
e  = np.array(e)
print(type(e))



u = [1,2,3,4]
u.reshape(-1,1)

class Nameless:
    def __init__(self, b, W, margin=1.e-6):
        """
        b should be a matrix whoese rows are the basis vectors of
        a vectorspace.

        W should be a matrix whose rows contain the apporeriate scalars
        such that when one of the rows is multipied with the matrix b,
        one should then get a matrix whose rwos, when summed, give you a
        vector that is a linear combination of the basis vectors.
        """

        if not isinstance(b, np.ndarray):
            raise TypeError("b should be a matrix")
        elif np.abs(ls.det(d)) < margin:
            raise Exception("basis vetcors are almost linearly dependent")

        if not isinstance(W, np.ndarray):
            raise TypeError("W should be a matrix")


        VectorMatrices =[]
        L = W.shape[1]
        for i in range(1, L):
            V = b *W[i]
            VectorMatrices.append(V)

        def Functional(b):
            if b != np.eye(L):
                raise Exception("b is not the identity matrix")
            else:
                Delta = b
                return Delta

        self.KroneckerDelta = Functional(b)

    def __CoordTransform__(self, X, XP, h=1.e-10):
        """
        This method is supposed to take the identity matrix of in tensor analysis
        mostly refered to as the "Kronecker delta" and transfrom it


        """

        def phi_psi(X,Xp):
            list_ =[]
            for x in X:
                x = x(*Xp)
                list_append(x)
            Epsilon = np.array(list_)
            return list_

        def fin_diff(f, x):
            """
            this function computes the derivative at a point x using finite difference.
            """
            fp = (f(x+h) -f(x))/h
            return fp

        def many_dot(reference, *M):
            """
            This function matrix-multplies the matrices in the list M

            PARAMETERS:
                reference (array):
                    this matrix is used to make sure that when we compare the elements
                    in M with it, the function raises an exception if one of the
                    elements in M are not the same dimension, by using M[0] as the reference.

                M (list):
                    this list contains the matrices that you want to multiply.

            RETURNS:
                Var[-1] (array):
                    Var[-1] is the resultant matrix.

            NOTICE:
                the order of the elements in M matter, matrix multiplication is genreally not
                commutative and this function runs in accordance with these rules.
            """
            if type(reference) != np.ndarray:
                raise ValueError(f"the reference, {reference} should be a an array")
            elif np.ndim(reference) != 2:
                raise ValueError("this function only handles matrices")
            else:
                Var =[np.eye(reference.shape[0])]
                for s in M:
                    if s.shape != reference.shape:
                        raise Exception(f"{M} contains a matrix whose shape differs from the reference")
                    new =np.matmul(Var[-1],s)
                    Var.append(new)
                return Var[-1]

        phi = phi_psi(X, Xp)
        L = len(phi)
        Jacobian = np.zeros((L,L))

        for x in phi:
            empty =[]
            for xp in Xp:
                x_diff = fin_diff(x, xp, h)
                empty.append(x_diff)
            r = np.array(empty)
            ind = phi.index(x)
            J[ind] = r


        M = [J.T, self.KroneckerDelta, J]
        G = many_dot(M[0], *M)
        self.G = G
        return Nameless(self.G)
























































class MultLinForm:
    def __init__(self,u0, *args):
        if not isinsatnce(u0, np.ndarray):
            raise TypeError("u0 should be a vector")
        """
        if isinstance(args, list):
            args = np.array(args)
            args = np.row_stack([u0, args])
            args = list(args)
        """

        for u in args:
            try:
                u0 *= u.reshape(-1,1)
            except AttributeError(f"{u} is a non-array type object"):
                u = np.array(u)
        self.args = args

        if len(self.args) == 0:
            LinearForm = u0
            self.G = LinearForm

        elif len(self.args) == 1:
            BilinearForm = u0
            self.G = BilinearForm
        else:
            MultinlinearForm = u0
            self.G = MultilinearForm

















        sigma = np.zeros_like(otherargs[0])
        for v in otherargs:
            if not isinstance(v, np.ndarray):
                v = np.array(v)
            sigma += W(v)*v


























"""
FIle I/O

    - working with measurded or scanned data
    - interacting with other programms
    - sawving information for comparisons or other postprosecessng needs
    FIle bject

    a file is a pythion obejct with associated metjods:
        # creating a read-only file object
"""
"""
myfile = open("measurement.data", "r")
"""
"""
   -  teh whole file can be read and stored into a string by:
"""
"""
s = myfile.read()
"""
"""

    - you can also read it one line at a time
"""
"""
myfile.readline()
"""
"""
    - or
"""
"""
for line in myfile.readlines():
    print(line)

"""
"""
FILES AS GENERATORS:


ITERABLE


the zip function couples elements in (any a ount of listsa actually) each lists with the same index into a tuple
list1 =[...]
list2 =[...]
list3 =[...]
zip(list1, list2, list3)
# returns a list with elements [(list1[0],list[0], list3[0]), (list[1], ...)... osv]

"""



"""
Making your own exception
you derive a class from a paren-class
"""
class FlatFunctionError(Exception):
    pass
    """
    the class FlatFunctionError , you don't need to add new methods
    to the derived class to be able to rasie completelly new Errors.
    WHat you do is tat when you are aware that you code will encounter
    an error that isn't built-in in pythin, you simply write
    """
class NUMA01Exception(Exception):
    pass

def f1(x):
    c = f2(x)
    return c
def f2(x):
    c = f3(x)
    return c
def f3(x):
    return 1/x

f3(0)






coord = [[1,1], [2,2], [2,2], [1,2], [0.5,1.5]]
coord.append(coord[0]) #repeat the first point to create a 'closed loop'

xs, ys = zip(*coord) #create lists of x and y values

plt.figure()
plt.plot(xs,ys)
plt.show() # if you need...
print(xs)
print(coord)


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



class PolyGon:
    def __init__(self,*points):
        l = len(points)
        if l < 3:
            raise Exception("you should provide self with at least three points")

        self.points = points
        self.points[0] = np.array(points[0])
        self.array = np.array(points)

        path = np.zeros_like(self.array)
        l = len(points)
        for i in range(1,l):
            path[i] = self.array[i]  - self.points[0]



        nonComplete_outer_hull = np.zeros((self.points[0].ndim, l+1))
        for i in range(1,l):
            nonComplete_outer_hull[i-1] = path[i] - path[i-1]



        outer_hull = np.row_stack([path[0], nonComplete_outer_hull, path[-1]])
        polyg = [list(s) for s in list(outer_hull)]
        PGarray=[s for s in list(outer_hull)]
        self.outer_hull = outer_hull
        self.polyg = polyg
        self.path =list(path)


    def area(self):
        def cross(A,B):
            if not isinstance(A, np.array) or isinstance(B, np.array):
                A = np.array(A)
                B = np.array(B)

            if A.ndim != 2:
                raise Exception("A should be in th eplane to for this class")
            elif A.ndim != B.ndim:
                raise Exception("A and B should both be in the plane")

            area = np.abs(A[0]*B[1] - A[1]*B[0])/2
            return area
        Area_list = []
        for s in self.path:
            if (s==self.path[-1]).all():
                return sum(Area_list)
            else:
                si = self.path.index(s)
                delta_area = cross(s,self.path[si +1])
                Area_list.append(delta_area)



    def plot(self):
        Coord = list(self.points)
        Coord.append(coord[0])
        xs, ys = zip(*coord)
        plt.figure()
        plt.plt(xs,ys)
        plt.show()

   # def contain(self, a):






 #   class Rectangle(PolyGon):

A = [(0.,0), (1.,0), (1.,1), (0.,1), (0.,0)]
AA = PolyGon(*A)














def integral(f,a,b,p, h=1.e-9, part=[], partp=[], partpp=[]):
    def n(p):
        return p
    def fp(x):
        diff_f = (f(x+h) - f(x))/h
        return diff_f


    part = [a + ((b-a)/n(p))*k for k in range(int(round(n(p)))+1)]
    partp = [fp(xi) for xi in part]
    partpp = [part[k]*partp[k] for k in range(len(part))]

    DataA = [f(xi)*(b-a)/n(p)/fp(xi) for xi in partpp]
    integral_value = sum(DataA)
    return integral_value



print(quad(np.sqrt, 0, 4))
print(integral(np.sqrt, 0,4, 100000))





def Time(p):
    t0 = time.time()
    def integral(f,a,b,p, h=1.e-9,N=None, part=[]):
        def n(p):
            return p

        part = [a + ((b-a)/n(p))*k for k in range(n(p) +1)]
        DataA = [f(xi)*(b-a)/n(p) for xi in part]
        integral_value = sum(DataA)
        return integral_value
    t1 = time.time()
    total = t1-t0
    return total

E = [Time(p) for p in range(1000)]
e = [p for p in range(1000)]
plt.plot(e, E, label ="execution time")
plt.xlabel=("p")
plt.ylabel=("time")
plt.title("comparing the error of the fast and slow approximation")
plt.legend(loc = "best")









t0 = time.time()
def scalarprod(A,B):
    if np.shape(A) != np.shape(B):
        raise Exception("Te two vectors should be of same dim")
    state = np.vstack([A,B])
    state2 = state.prod(axis=0)
    scalar = state2.sum(axis=0)
    return scalar
M = np.array([[1,1,1,1],
             [2,2,2,2]])

print(scalarprod(M[0], M[1]))
t1 = time.time()
total = t1-t0
print(total)

t0 = time.time()
def vv(A, B):
    if type(A) != type(B) != np.ndarray:
        A = np.array(A)
        B = np.array(B)
    elif len(A) != len(B):
        raise Exception("your input vectors should have the same dimesion")
    return sum(list(A * B))
M = np.array([[1,1,1,1],
             [2,2,2,2]])
print(vv(M[0], M[1]))
t1 = time.time()
total = t1-t0
print(total)

M[1]



def scalarprod(A,B):
    if np.shape(A) != np.shape(B):
        raise Exception("Te two vectors should be of same dim")
    state = np.vstack([A,B])
    state2 = state.prod(axis=0)
    scalar = state2.sum(axis=0)
    return scalar







"""
lecture 7
reshape:
    takes a vector and makes of it a matrix.
    A = arange(6) # arrqy with elements from 0-5
    A.rehspae(6,1)
    A.reshpape(2,3)
    A.reshape(3,2)
    A.rehsape(1,6)

    Rehsape has the clever way of working like "i give you -1 and python gives you
    a shepa that auotomatically fits"
    Exampe


    a = arange(12)
    A.reshape(1,-1) # row matrix
    A .reshape(-1,1) # cloumn matrix
    A. reshape(3,-1) # shape = 3,4
    A. rehsape(-4,-1)

STacking Vetcors:
    v1 = np.array([1,2])
    v2 = np.array([3,4])

    M = vstack([v1,v2])
    N = column_stack([v1,v2])

SUM, MAX, MIN:
    -  A.sum()
        if i don't give a parameter in the bracket, the fuctions sums every
        elemnt in the array.
    - A.sum(axis=0)
    - A.min




"""
M = np.array([[1,1,1,1],
             [2,2,2,2]])
print(M.sum(axis=0))
print(M.sum(axis=1))
print(M.prod(axis=0))


"""

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



def fastlogapprox_amer(x,n, fast_diff_mat="off"):

    a,g = [(1+x)/2], [np.sqrt(x)]

    for i in range(n):
        a_new = (a[-1] + g[-1])/2
        a.append(a_new)
        g_new = np.sqrt(g[-1]*a[-1])
        g.append(g_new)

    a = np.array(a)
    d[:1] = a
    for j in range(1,n+1):
        for i in range(j,n+1):
            d[j,i] = (d[j-1,i] -(d[j-1,i-1]*4**(-j)))/((1-4**(-j)))

    if fast_diff_mat == "on":
        d[-1:,0] = np.abs((x-1)/d[-1,-1] - (x-1)/a[-1])
        return d
    elif fast_diff_mat=="off":
        return d,(x-1)/d[-1,-1]


print(fastlogapprox_hector(5, 5))
print("-------------------------------------------------------------")
print(fastlogapprox_alex(5, 5))
print("-------------------------------------------------------------")
print(fastlogapprox_amer(5, 5))


"""



"""
Lecture 6?

PARTIAL APPLICATION:
    making a fucntion that returns another function:

        def make_sine(omega):
            def f(x):
                return np.sin(omega*x)
            return f
        - f is a function of a single parameter, x
        - make_sine creates a sine function for every value of omega.

    There is another way of using partial application: (Anonoymous functions)

        -example:
            f = lambda x: x**3 + x*2 +1
            g = lambda x,y : 3*x-2*y
            quad(lambda x: <någon_function>(x), a,b)
        - make sine function can now be written as:
            omega =3
            fomega = lambdax: np.sin(omega*x)
            si.quad.()

        SÖK FÖR import functools

    ARRAYS:
        arrays take only components of the same datatype. if you'd forget the point which
        makes an integer turn into a float, then numpy makes all the elements into a float,
        aka, it makes the whole array into a "same data type"-array

        example

        A = np.array([1,2., 3.,]) turns into np.array([1.,2.,3.])
        B = np.array([1.,2.,1+ 2*1*j])

        1/2    #0.5
        1//2   # rounds donw to nearest integer , 0
        3%2   # is the modulus´, gives the rest, 1
        - creating matrices:
            - full((n,m),a) returns a ,matrix of shape n,m with the entries




"""
a = 1+ 2*1j
print(type(a))



b = np.empty(4)
print(b)
















"""
def R(x,y):
    return np.sin(x)*np.cos(y)
def I(x,y):
    return x*y
"""

class Complex:
    def __init__(self, R,I, x,y, re=None, im=None):
        if re == 0:
            self.re = lambda x,y: 0
        elif re == x:
            self.re = lambda x,y: x
        elif not re:
            self.re = R(x,y)

        if im == 0:
            self.im = lambda x,y: 0
        elif im == y:
            self.im = lambda x,y: y
        elif not im:
            self.im = I(x,y)




    def __repr__(self):
        if self.im < 0:
            return f"{self.re} - {(-1)*self.im}i"
        else:
            return f"{self.re} + {self.im}i"
    """-------------------------------------------------------------------------------------
    """


    def __add__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.re + other[0].re, self.im + other[0].im)
        elif len(other) > 1:
            val1, val2 = self.re, self.im
            for s in other:
                val1 += s.re
                val2 += s.im
        return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """
    def __sub__(self, *other):
        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return Complex(self.re - other[0].re, self.im - other[0].im)
        elif len(other) > 1:

            val1, val2 = self.re, self.im
            for s in other:
                val1 -= s.re
                val2 -= s.im
            return Complex(val1, val2)

    """------------------------------------------------------------------------------------
    """
    def __mul__(self, *other):
        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)

        if len(other) == 0 :
            return 0
        elif len(other) == 1:
            return Multiply(self, other[0])
        elif len(other) > 1:
            start_self = [self]
            for s in other:
                new = Multiply(start_self[-1], s)
                start_self.append(new)
                return start_self[-1]

    """------------------------------------------------------------------------------------
    """
    def __div__(self, *other, margin = 1.e-10):
        def Multiply(a,b):
            """
            this function multiplies two instances of the class Complex
            """
            a1, a2 = a.re, a.im
            b1, b2 = b.re, a.im
            return Complex(a1*b1 -a2*b2, a1*b1 + a1*b1 +a2*b2)
        def conj(a):
            a1, a2 = a.re, a.im
            return Complex(a1, -a2)

        def Divide(a,b):
            """
            this function divides thw instances of the class Complex
            """

            b_conj = conj(b)
            c = Multiply(a, b_conj)
            absb = float(b_conj.re)**2 + float(b_conj.im)**2
            return Complex(c.re/absb, c.im/absb)

        start_self = [self]
        for s in other:
            abss = float(s.re)**2 + float(s.im)**2
            if abss < margin:
                raise Exception("div, other = {other} contains a number that is too small")
            new = Divide(start_self[-1], s)
            start_self.append(new)
        return start_self[-1]

    """---------------------------------------------------------------------------------------
    """

    def __get_attr__(self, *other):
        def prel_attr(self,*other, Re=[], Im=[]):
            Re = [s.re for s in other]
            Im = [s.im for s in other]
            return [Re, Im]

        if len(other) == 0 :
            return self
        elif len(other) == 1:
            return prel_attr(self, other[0])
        elif len(other) > 1:
            return prel_attr(self,*other)

    """------------------------------------------------------------------------------------------
    """
    def __is_attr__(self, *other):
        """
        I do't think this method is very good. Needs revision.
        """
        def ineq_(a, c):
            def is_real_prel(a):
                if a.im == 0:
                    return True
                else:
                    return None
            def is_imag_prel(a):
                if a.re == 0:
                    return True
                else:
                    return None
            if np.abs(float(a.re)) < np.abs(float(a.im)):
                c ==1
                if c==1:
                    return is_real_prel(a,c)
            elif np.abs(float(a.re)) > np.abs(float(a.im)):
                c == 2
                if c ==2:
                    return is_imag_prel(a,c)

        empty = []
        for s in other:
            empty.append(ineq_(s,c))
        return empty

def f(x,y):
    return np.cos(x*y)
def g(x,y):
    return np.sin(x*y)

x,y =2,7
A = Complex(lambda x,y: np.cos(x), lambda x,y: np.sin(y),  x,y)
print(A)
AA = [A]
a = Complex.__get_attr__(AA)
print(a)
X =np.linspace(0,np.pi/2,100)
Q = [Complex(lambda i,y: np.cos(i), lambda x,i: np.sin(i), ) ]
print(A.re)

"""
S = Complex(3,3)
Q = Complex(1,2)
P = Complex(2,3)
other = [Q, P, S]
a = Complex.__add__(*other)
b =Complex.__sub__(*other)
c =Complex.__mul__(*other)
d = Complex.__div__(*other)
e = Complex.__get_attr__(*other)
f = Complex.__is_attr__(*other)
print(I)
print([a,b,c,d,e,f])

"""

class Interval:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __repr__(self):
        return f"[{self.a}, {self.b}]"

    def __add__(self, *other):
        v1, v2 = self.a, self.b
        if len(other)==0:
            return self
        elif len(other)==1:
            return Interval(self.a +other[0].a, self.b + other[0].b)
        elif len(other) > 1:
            v1, v2 = self.a, self.b
            for s in other:
                v1 += s.a
                v2 += s.b
            return Interval(v1, v2)

    def __sub__(self, *other):
        v1, v2 = self.a, self.b
        if len(other)==0:
            return self
        elif len(other)==1:
            return Interval(self.a - other[0].b, self.b - other[0].a)
        elif len(other) > 1:
            v1, v2 = self.a, self.b
            for s in other:
                v1 -= s.b
                v2 -= s.a
            return Interval(v1, v2)

    def __mult__(self, *other):
        def Multiply(A,B):











































