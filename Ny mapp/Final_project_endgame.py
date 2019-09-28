# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 13:55:02 2019

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

E = np.array([1.,2])
M = np.array([[1.,1],
             [4,5]])
print(np.matmul(E, M))


a = (1/1.6605/10**(-27))*1.602*10**(-19)/1000000
print(a**2)

def delta_M(A, Z, atomic_mass , EB=True, eV=True, me=9.10938356*10**(-31), mp=1.67262192369*10**(-27), mn=1.674927*10**(-27), c=3*10**(8), u=1.6605*10**(-27), e=1.602*10**(-19)):
    mp=1.67262192369*10**(-27)/u
    mn=1.674927*10**(-27)/u
    me=9.10938356*10**(-31)/u
    N = A-Z

    MH= 1.007825
    dM = Z*MH + N*mn - atomic_mass
    return dM, mp, mn, me, atomic_mass

def energy(args, argsp,c=3*10**(8) ):
    Q = sum(args) - sum(argsp)*c**2
    return Q



dm = delta_M(62, 28, 61.928345)









def E_B(A, Z, eV=True, me=9.10938356*10**(-31), mp=1.67262192369*10**(-27), mn=1.674927*10**(-27), c=3*10**(8), u=1.6605*10**(-27), e=1.602*10**(-19)):
    mp=1.67262192369*10**(-27)/u
    mp=1.67262192369*10**(-27)/u
    me=9.10938356*10**(-31)/u


    N = A-Z
    def mass(A,Z):
        N = A-Z
        return Z*(me+mp) + N*mn

    MH= 1.007825 + me
    EB = (Z*MH + N*mn - mass(A,Z))*c**2
    if eV:
        return EB/e
    else:
        return EB


def mass_deffect(A, Z, c=3*10**(8)):
    return E_B(A, Z)/c**2

print(mass_deffect(62, 28))






def E(n):
    h=6.63*10**(-34)
    e=1.602*10**(-19)
    c=3*10**(8)
    nano=1*10**(-9)
    me=9.11*10**(-31)
    eps = 8.854*10**(-12)
    K = me*e**4/(8*eps**2*h**2)
    return 1/(-K*n**2)

fig, ax = plt.subplots()
N = np.array(range(2, 101))
N_ = np.array(range(1, 100))
NN = np.array([[N[-i] -N_[j] for i in range(N.shape[0])] for j in range(N_.shape[0])])
for k in range(NN.shape[0]):
    ax.plot(NN.T[k], E(NN.T[k]), label = "energy levels of hydrogen atom")


fig, ax = plt.subplots()
N = np.array(range(1, 100))
ax.plot(N, E(N), label = "energy levels of hydrogen atom")


a = e - h*c/600
print(a)

def direct_projection(self, T1, T2, p_a1, p_a2, contraction_order=1, p_a3=None):
    if contraction_order == 3 and not p_a3:
        raise Exception("you must specify p_a3 to be able to perform a self projection on a rank 4 tensor")
    else:
        def m(a,b, out=None):
            return np.matmul(a,b, out=None)
        def Gtr(a,b, M=None, out=None, M_direct=True):
            """
            M is an already formed matrix
            """
            if not M:
                return np.trace(m(m(a,b), G))
            elif not a and not b:
                return np.trace(m(M, G))
            elif direct:
                return np.trace(M)

        r1 = np.ndim(T1)
        r2 = np.ndim(T2)
        def TG_primer(T, p_a, rank):
                """
                THis function primes the tensor for the lowering or the riase
                of an index by the the function TG bellow. aka, it chooses
                which index TG performs index juggling on.
                """
                def rotation_module(T, p_a, rank=rank):
                    A = T
                    if rank == 3:
                        d_c = [[-1,0,0],
                               [0,-1,0],
                               [0,0,-1]]
                        L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]
                        if p_a == d_c[0]: #the motivation for this is that: imagine ghat you had an axis in the negative x-direction attached to the cube
                                          # and thta you rotate clockwise, the face facing the forward direction will now
                                          #be facing the uppward direction, that is why the non-zero entries are -1.
                                          # the same logic applies for the next elif command but in the negative y-axis deriction.
                                          # the actuall programming process has nothing to do with this visualisation
                            A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                            return A_T
                        elif p_a == d_c[1]: # refoemration parallelel to y-axis
                            A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                            return A_T
                        elif p_a == d_c[2]: #leaves it be
                            return A
                    elif rank == 2:
                        d_c = [[-1,0],
                               [0,-1]]
                        if p_a == d_c[0]:
                            return A.T
                        elif p_a == d_c[1]:
                            return A
                    elif rrank == 1:
                        return A

                e =  [1,2,3]
                if rank in e:
                    return rotation_module(T, p_a)
                elif rank == 4:
                    d_c  = [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]
                    def rank4_rotation(T, t):
                        def upward_shelve(T):
                            L = np.array([[T[i, j] for i in range(T.shape[0])] for j in range(T.shape[1])])
                            return L
                        Tr = np.array([rotation_module(T[i], t) for i in range(T.shape[0])])
                        L = upward_shelve(Tr)
                        return L
                    if p_a == d_c[0]:
                        return rank4_rotation(T, [-1,0,0])
                    elif p_a == d_c[1]:
                        return rank4_rotation(T, [0,-1,0])
                    elif p_a == d_c[2]:
                        return rank4_rotation(T, [0,0,-1])
                    elif p_a == d_c[3]:
                        return T

        def TT_naive(T1, T2, p_a1, p_a2, r1, r2, contraction_order):
            T1 = TG_primer(T1, p_a1, r1)
            T2 = TG_primer(T2, p_a2, r2)

            e = [1,2]
            def MV_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order):
                if contraction_order == 1:
                    T1 = TG_primer(T1, p_a, r1)
                    return m(T1, T2)
                elif r1 != 2 and r2!= 2 and contrection_order == 2:
                    raise Exception("you cant contract twice if their first order contraction is a vector")
                elif contrection_order == 2:
                    T1 = TG_primer(T1, p_a1)
                    T2 = TG_primer(T2, p_a2)
                    T3 =  m(T1, T2)
                    s  =  Gtr(a=None, b=None, M=T3)
                    return s

            def TT_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order, p_a3=p_a3):
                def order4_trace(T1, T2, p_a3=p_a3):
                    L = np.array([[m(T1[i], T2[j]) for i in range(T1.shape[0])] for j in range(T2.shape[0])])
                    r3 = np.ndim(L)
                    L = TG_primer(L, p_a3, r3)
                    M = np.array([[Gtr(a=None, b=None, M=L[i,j]) for i in range(L.shape[0])] for j in range(L.shape[1])])
                    return M

                if contraction_order==1:
                    if r2 == 1:
                        T1 = TG_primer(T1, p_a1, r1)
                        L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                        return L
                    else:
                        T2 = TG_primer(T2, p_a2, r2)
                        if r2 ==2:
                            L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                            return L
                        elif r2==3:
                            L = np.array([[m(T1[i], T2[j]) for i in range(T1.shape[0])] for j in range(T2.shape[0])])
                            return L
                elif contraction_order == 2:
                    if r2 ==1:
                        T1 = TG_primer(T1, p_a1, r1)
                        L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                        s = Gtr(a=None, b=None, M=L)
                        return s
                    elif r2 == 3:
                        return order4_trace(T1, T2)

                elif contraction_order == 3:
                    a = order4_trace(T1, T2)
                    s = s = Gtr(a=None, b=None, M=a)


            if r1 in e and r1 in e:
                return MV_comb(T1, T2, p_a1, p_a2, rank1, rank2)
            else:
                return TT_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order)

        return TT_naive(T1, T2, p_a1, p_a2, r1, r2, contraction_order)























print(np.sqrt(2*(1.602*10**(-19)*1.25)/9.11/10**(-31)))








def TG(self,T, axis, L_R):
        """
        this function lowers or raises the index of a tensor of rank 1,2, or 3
        axis specifies which index and is a list
        L_R stands for Lower_Raise
        """
        rank = np.ndim(A)
        ttt = Tensor_type(self, T, args)
        p_a = axis
        def m(a,b, out=None):
            return np.matmul(a,b, out=None)

        def TG_primer(T, p_a, rank=rank):
            """
            THis function primes the tensor for the lowering or the riase
            of an index by the the function TG bellow. aka, it chooses
            which index TG performs index juggling on.
            """
            def rotation_module(T, p_a, rank=rank):
                A = T
                if rank == 3:
                    d_c = [[-1,0,0],
                           [0,-1,0],
                           [0,0,-1]]
                    L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]
                    if p_a == d_c[0]: #the motivation for this is that: imagine ghat you had an axis in the negative x-direction attached to the cube
                                      # and thta you rotate clockwise, the face facing the forward direction will now
                                      #be facing the uppward direction, that is why the non-zero entries are -1.
                                      # the same logic applies for the next elif command but in the negative y-axis deriction.
                                      # the actuall programming process has nothing to do with this visualisation
                        A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                        return A_T
                    elif p_a == d_c[1]: # refoemration parallelel to y-axis
                        A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                        return A_T
                    elif p_a == d_c[2]: #leaves it be
                        return A
                elif rank == 2:
                    d_c = [[-1,0],
                           [0,-1]]
                    if p_a == d_c[0]:
                        return A.T
                    elif p_a == d_c[1]:
                        return A
                elif rank == 1:
                    return A

            e =  [1,2,3]
            if rank in e:
                return rotation_module(T, p_a)

            elif rank == 4:
                d_c  = [[-1,0,0,0],
                        [0,-1,0,0],
                        [0,0,-1,0],
                        [0,0,0,-1]]

                def rank4_rotation(T, e):
                    def upward_shelve(T):
                        L = np.array([[T[i, j] for i in range(T.shape[0])] for j in range(T.shape[1])])
                        return L
                    Tr = np.array([rotation_module(T[i], e) for i in range(T.shape[0])])
                    L = upward_shelve(Tr)
                    return L

                if p_a == d_c[0]:
                    return rank4_rotation(T, [-1,0,0])
                elif p_a == d_c[1]:
                    return rank4_rotation(T, [0,-1,0])
                elif p_a == d_c[2]:
                    return rank4_rotation(T, [0,0,-1])
                elif p_a == d_c[3]:
                    return T

        def Lower_Raise(T, p_a, L_R, ttt=ttt):
            def G_G(self, L_R):
                if L_R == [-1,1]:
                    return self.metric_tensor
                elif L_R == [1,-1]:
                    return self.inv_metric_tensor

            if L_R == [-1,1] or L_R == [1,-1]:  #----> lowering & raising indices respectivelly
                if ttt[1] == 0 and L_R[1] == -1 or ttt[0] == 0 and L_R[0] == -1:
                    raise Exception(f"you cant raise/lower an already raised/lowered tensor index. Your tensortype {ttt}")
                else:
                    T = TG_primer(T, p_a)   #------> this is the part where p_a is used

                    if rank == 1:
                        L = m(G_G(L_R), T)
                    elif rank == 2:
                        L = np.array([m(G_G(L_R), T)])
                    elif rank == 3:
                        L = np.array([m(G_G(L_R), T[i]) for i in range(T.shape[0])])

                    ttt = [ttt[0] +L_R[0], ttt[1] +L_R[1]]
                    return L,ttt
            else:
                raise Excperion("the following must hold: Lower_raise == [-1,1] or Lower_Riase == [1,-1]")
        return Lower_Raise(T, p_a, L_R)


def direction_catalogue(self,rank):
    d_c4  = [[-1,0,0,0],
            [0,-1,0,0],
            [0,0,-1,0],
            [0,0,0,-1]]

    d_c3  = [[-1,0,0,0],
             [0,-1,0,0],
             [0,0,-1,0],
             [0,0,0,-1]]

    d_c2 = [[-1,0],
            [0,-1]]
    if rank ==4:
        return d_c4
    elif rank ==3:
        return d_c3
    elif rank ==2:
        return d_c2

























































a =1
b = 2
print(not a and b)


H = np.array([[0., 1.,  2.,  33.],
              [1., 5.,  6.,  66.],
              [2., 9.,  10., 84.],
              [3., 12., 13., 99.],
              [4., 15., 16., 46.],
              [5., 18., 19., 87.],
              [6., 21., 22., 45.],
              [7., 24., 25., 77.],
              [8., 27., 28., 44.]])

print(H.T[0])
a = np.array([1,3,2,5])
b = np.array([2,4,2,7])
c = b - a.reshape(-1,1) +np.eye(4)
e = np.array([3,9,4,1])
d = np.array([e.reshape(-1,1)*c[i] for i in range(4)])
print(d)

def Tensor_index_rotator(A, axis, lin_comb_axis=None, rank=None):
    """
    LOL,This function is so good that you can just input a tensor twice and ot will revert back
    to its original state, just like a matrix transpose fucntion.
    direction_list can only be [-1,0,0], [0,-1,0] or [0,0,-1] for rank=3,
    elsedirection_list =NOne
    To do the reverse transpose of the transpose you've made, just turn
    the default value for inverse to True
    """
    direction_catalogue = [[-1,0,0],
                           [0,-1,0],
                           [0,0,-1]]
    L_A1,L_A2,L_A3 = A.shape[0], A.shape[1], A.shape[2]

    if np.ndim(A) == 2:
        return A.T

    elif np.ndim(A) ==3:
         def transpose_machinery(A, direction_list):
             if axis == direction_catalogue[-1]:
                 A_T = np.array([A[i].T for i in range(L_A1)])

             elif axis == direction_catalogue[-2]:
                 A_T = np.array([np.array([A[j].T[i] for j in range(L_A1)]).T for i in range(L_A3)])

             elif axis == direction_catalogue[-3]:
                 A_T = np.array([np.array([A[j,i] for j in range(L_A1)]).T for i in range(L_A2)])

             return A_T

         if not lin_comb_axis:
             return transpose_machinery(A, direction_list)
         elif lin_comb_axis:
             if L_A1 == L_A2 == L_A3:
                 l = len(direction_catalogue)
                 L = [transpose_machinery(A, direction_catalogue[k]) for k in range(l)]
                 new = L[0] + L[1] + L[2]
                 return new
             else:
                 raise Exception(f"your tensor must be a perfect cube,yours has the dime: {K}")
C = Tensor_index_rotator(d, [0,0,-1])
print("------------------LOOL1-----------------------")
print("------------------LOOL2-----------------------")
print(Tensor_index_rotator(C, [0,0,-1]))








#-------------------------------------------------------------------------------------------------------------------------------------------



def contravariant_MM_projector(self, M1, M2, T1=None, T2=None, trace=False):
    """
    Ti is a reference to the mathematical operation of taking one of the indices of the
    matrix Mi and contracting it wit an index of the metric tensor. if the index of Mi is the
    first index then them matrix Mi is not transoosed and Ti = None. Else, if one wishes
    to use the other index as means of projecting th metric tensor on the matrix.
    this is basically the same as lowering one index, either the first or the second index.
    """
    G = self.metric_tensor
    G_ = self.inv_metric_tensor

    def m(a,b, out=None):
        return np.matmul(a,b, out=None)
    def Gtr(a,b, M=None, out=None):
        """
        M is an already formed matrix
        """
        if not M:
            return np.trace(m(m(a,b), G))
        elif not a and not b:
            return np.trace(m(M, G))

    def Conditional_transpose(M1,M2):
        E =[T1, T2]
        if not T1 and not T2:
            return (M1, M2, (0,0))
        elif not T1 or not T2:
            ind = E.index(None)
            if ind ==0:
                return (M1, M2, (0,1))
            else:
                return (M1, M2, (1,0))

    def T_verif(a, indication_number):
        n = indication_number
        if n == 0:
            return a
        elif n == 1:
            return a.T

    def final_step(M1, M2):
        a,b,indicator = *Conditional_transpose(M1,M2)
        if indicator == (0,0):
            return m(a, m(G,b))
        else:
            obj = m(T_verif(a,indicator[0], m(G, T_verif(b, indicator[1]))))
            return obj

    if trace == False:
        return final_step(M1,M2)
    else:
        Gtr(None, None, final_step(M1,M2))



#------------------------------------------------------------------------------------------------------------------------------------


def co_contra_MM_projector(self, M1, M2, p_a1, p_a2, L_R1, L_R2):
    ttM1, ttM2 = Tensor_type(self, M1, args), Tensor_type(self, M2, args)
    d_c = [[-1,0], [0,-1]] # this list codes for the first index and second index repsectivelly of a rank 2 tensor
    G = self.metric_tensor
    G_ = self.inv_metric_tensor

    def m(a,b, out=None):
       return np.matmul(a,b, out=None)


    def TG(T, p_a, ttt, L_R):
        """
        this function lowers or raises the index of a tensor
        """
        def TG_primer(a, p_a):
            """
            THis function primes the tensor for the lowering or the riase
            of an index by the the function TG bellow. aka, it chooses
            which index TG performs index juggling on.
            """
            if rank == 3:
                L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]

                if p_a == d_c[0]: #reformation parallel to x-axis
                    A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                    return A_T

                elif p_a == d_c[1]: # refoemration parallelel to y-axis
                    A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                    return A_T

                elif p_a == d_c[2]: #leaves it be
                    return A

            elif rank == 2:
                if p_a == d_c[1]:
                    return A
                elif p_a == d_c[0]:
                    return A.T

        if L_R == [-1,1] or L_R == [1,-1]:
            if ttt[1] == 0 and L_R[1] == -1 or ttt[0] == 0 and L_R[0] == -1:
                raise Exception(f"you cant raise/lower an already raised/lowered tensor index. Your tensortype {ttt}")
            else:
                def Low_Rise(T, p_a, L_R):
                    if L_R == [-1,1]:
                        G = G
                    elif L_R == [1,-1]:
                        G = G_
                    T = TG_primer(T)
                    L = np.array([MV_(T[i],G[i]) for i in range(T.shape[0])])
                    ttt = [ttt[0] +L_R[0], ttt[1] +L_R[1]]
                    return L,ttt
                return Low_Rise(T, p_a, L_R)[0]
        else:
            raise Excperion("the following must hold: Lower_raise == [-1,1] or Lower_Riase == [1,-1]")

    M1, M2 = TG(M1, p_a1, ttM1, L_R1), TG(M2, p_a2, ttM2, L_R2)
    return m(M1, M2)



#-------------------------------------------------------------------------------------------------------------------------------------------








def co_contra_TV_projector(self,T, v, primer_axis, Lower_Raise):
    """
    this function (hopefully without any ptoblems) projects a vector onot a rank 3 tensor.

    inputs:
        T (tensor)
        V (vector)
        primer_axis (list):
            you can choose between the sub-lists inside the matrix shpaed list
            bellow (direction_catalogue)
        Lower_Raise (list):
            either [-1,1] or [1,-1] for lowering and raising indices respectivelly

    Returns:
        L (tensor):
            either a vector or a mtrix. I should expand this a bit further.

    The problem with this fucntion is that i'm contructing it in a somewhat hazy state of
    thinking. I'm beginning to feel like this is going to become a bit obsessive, just like
    when i was obsessing over einsteins equations and trying to derive his field equations.
    Fuck.

    """
    p_a     = primer_axis
    ttt     = Tensor_type(self, T, args)
    ttv     = Tensor_type(self, v, args)
    G       = self.metric_tensor
    G_      = self.inv_metric_tensor
    rank    = np.ndim(T)
    if rank == 3:
        direction_catalogue = [[-1,0,0],
                               [0,-1,0],
                               [0,0,-1]]
    elif rank == 2:
        direction_catalogue = [[-1,0],
                               [0,-1]]
    d_c = direction_catalogue

    def TV_naive(T, v, p_a):
        T = TG_primer(T, p_a)
        L = np.array([m(T[i], v) for i in range(T.shape[0])])
        return L

    def m(a,b, out=None):
        return np.matmul(a,b, out=None)

    def TG_primer(A, p_a=primer_axis, rank=rank):
        """
        THis function primes the tensor for the lowering or the riase
        of an index by the the function TG bellow. aka, it chooses
        which index Tg performs index juggling on.
        """
        if rank == 3:
            L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]

            if p_a == d_c[0]: #reformation parallel to x-axis
                A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                return A_T

            elif p_a == d_c[1]: # refoemration parallelel to y-axis
                A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                return A_T

            elif p_a == d_c[2]: #leaves it be
                return A
        elif rank == 2:
            if p_a == d_c[1]:
                return M
            elif p_a == d_c[0]:
                return M.T


    def TG(General_Tensor, p_a=primer_axis, ttt=ttt, Lower_Raise=Lower_Raise):
        """
        this function lowers or raises the index of a tensor
        """
        T = General_Tensor
        L_R = Lower_Raise
        if L_R == [-1,1] or L_R == [1,-1]:
            if ttt[1] == 0 and L_R[1] == -1 or ttt[0] == 0 and L_R[0] == -1:
                raise Exception(f"you cant raise/lower an already raised/lowered tensor index. Your tensortype {ttt}")
            else:
                def Low_Rise(T, p_a, L_R):
                    if L_R == [-1,1]:
                        G = G
                    elif L_R == [1,-1]:
                        G = G_
                    T = TG_primer(T)
                    L = np.array([MV_(T[i],G[i]) for i in range(T.shape[0])])
                    ttt = [ttt[0] +L_R[0], ttt[1] +L_R[1]]
                    return L,ttt
                return Low_Rise(T, p_a, L_R)[0]
        else:
            raise Excperion("the following must hold: Lower_raise == [-1,1] or Lower_Riase == [1,-1]")

    T = TG(T)
    L = TV_naive(T, v)
    return L

#-------------------------------------------------------------------------------------------------------------------------------------------


















































E = np.array([[1,2,3],
              [4,5,6],
              [7,8,9]])

t = np.array([0,0,0])
C = np.vstack([E, t])
print(C.shape)
a = 1,0
b = 2,3
L = [(a[i], b[j]) for i in range(2) for j in range(2)]
print(L)


#D = np.array([e.reshape(-1,1)*d[i] for i in range(3)])
print(c)
print("-------------------------")
print(d)
print("-------------------------")
print("-------------------------")
print(D)

#x, y = symbols('x y', real=True)
#f = y*(x/(1+x)) +1 +x
#F1 = 3*x**3-3*x*y**2 -1
#F2 = 3*x**2*y -y**3
#
#a,b= diff(F1,x), diff(F1,y)
#c,d= diff(F2, x), diff(F2, y)

def n(x,y):
    return 3*x**3-3*x*y**2 -1

def deriv(x,y,n):
    x, y = symbols('x y', real=True)
    return diff(n,x)

E = np.array([[1,2,3],
              [4,5,6],
              [7,8,9]])
t = np.array([0,0,0])
C = np.vstack([E, t])
print(C.shape)
print(np.matmul(C, C.T))











class Manifold:
    def __init__(self, X, J, H=None):
        """
        You need to initialize this class with parameterised representation of a manifold.
        think of spherical coordinates that plot out points on a spherical shell for each value of
        the radius. each shell is its own manifold
        """
        if not isinstance(X, np.ndarray) and isinstance(X, list) or isinstance(X, tuple):
            X = np.array(X)
        else:
            raise TypeError("X should be a vector or be convertable to a vector")

        if not isinstance(J, np.ndarray) and isinstance(J, list) or isinstance(J, tuple):
            J = np.array(J)
        else:
            raise TypeError("X should be a vector or be convertable to a vector")
        self.Hermitian = H
        self.X = X
        self.B = B
        self.J = J
        Jt = J.T
        self.Jt = Jt
        G = np.matmul(J.T, J)
        self.metric_tensor = G
        self.inv_metric_tensor = np.linalg.inv(G)


    def Tensor_rank(self, T):
        rank = np.ndim(T)
        return rank


    def univ_finite_diff(self, Tensor, args, h=1.e-9):
        T = Tensor
        rank = np.ndim(T)

        def normal_deriv(T):
            """
            This function calculates the usual partial derivative of a tensor valued function.
            It is divided into nested modules that, depending on the rank of the tensor-input,
            call deeper modules to do the job. The upper-most module turns vectors into jacobians of sorts,
            the second module uses the first module to turn the rows of a rank 2 tensor into a 3D-cube of derivates,
            and so on.

            INPUT:
                T (np.ndarray):
                    tensor of maximum rank 3
            RETURNS:
                tensor of one rank higher than T containing information about how it changes from point to point.
                the output tensor has a rank equal to rank(T)+1
            """
            def deriv_module1(A):
                """
                A must be a vector
                """
                L_a     = A.shape[0]
                vargs   = 2*np.array(args)
                n       = vargs.shape[0]
                M       = vargs-0.5*vargs.reshape(-1,1) + h*np.eye(n)
                L_m     = M.shape[0]

                Jacobi = np.array([[(A[i](M[j]) - A[i](args))/h  for j in range(L_m)] for i in range(L_a)])
                return Jacobi

            def deriv_module2(B):
                L_b = B.shape[0]
                Gamma = np.array([deriv_module1(T[i]) for i in range(L_b)])
                return Gamma

            def deriv_module3(C):
                L_c = C.shape[0]
                Riemann = np.array(deriv_module2(T) for i in range(L_c))
                return Riemann

            if rank == 1:
                return deriv_module1(T)
            elif rank == 2:
                return deriv_module2(T)
            elif rank == 3:
                return deriv_module3(T)
            else:
                raise Exception(""""this class has no code for handling derivatives of
                                tensors of rank greater than 3, but it is easily fixable""")
        return normal_deriv(T)

    def surface2D_curvature_tensor(self, args):
        """
        This method calculates the curvature tensor for a surface using local
        gaussian normal coordinates and the normal itself. This method is only works
        on 2D surfaces, this is due to the fact that the cross product only works when
        vectors live an Ambient space of maximum 3 dimensions.

        this tensor {b_myu_ny} is relatied to the riemann curvature tensor
        in a special way. I might add that as a methdod too.

        INPUT:
            args (list, vector, tuple):
                i wasn't really sure how i would implement the finite difference
                without a coordinate system for the inputs to the tensor object.
                one needs coordinates to be able to differentiate function. So, the args input
                should contain coordinates for a point on the surface (manifold) using curvilinear
                coordinates.
        RETURN:
            b (matrix):
                containing
        """
        J = self.J

        if J.shape[1] == 2:
            def cross_prod(a,b, normalised=True):
                A = np.array([[0, -a[2], a[1]],
                              [a[2], 0, -a[0]],
                              [-a[1], a[0], 0]])
                C = np.matmul(A, b)
                if normalised:
                    return C
                else:
                    c = C/np.linalg.norm(C)
                    return c
            n = cross_prod(J[0], J[1])
            n_ny = self.univ_fifnite_diff(n, args)

            curvature_tensor = np.matmul(J.T, n_ny)
            b = curvature_tensor

            return b
        else:
            raise Exception("this tensor works ony on 2D surfaces")


    def Newt_motor(self, v0, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, finite=False, once=False):
            F       = self.X
            J       = self.J
            xvals   = [v0]
            x_max   = 1/margin
            k       = 0

            def solving_motor(v):
                def univ_finite_diff(v, h):
                    self.J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                                       [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                    return self.J

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
                        delta = np.linalg.solve(univ_finite_diffe(v, h), -F(v))
                        new = delta + v
                        return new
                    elif once == True:
                        J_once = univ_finite_diff(v0, h)
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

    def TG_primer(self,T, p_a, rank):
        """
        THis function primes the tensor for the lowering or the riase
        of an index by the the function TG bellow. aka, it chooses
        which index TG performs index juggling on.
        """
        def rotation_module(T, p_a, rank=rank):
            A = T
            rrank = np.ndim(A)
            if rrank == 3:
                d_c = [[-1,0,0],
                       [0,-1,0],
                       [0,0,-1]]
                L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]
                if p_a == d_c[0]: #the motivation for this is that: imagine ghat you had an axis in the negative x-direction attached to the cube
                                  # and thta you rotate clockwise, the face facing the forward direction will now
                                  #be facing the uppward direction, that is why the non-zero entries are -1.
                                  # the same logic applies for the next elif command but in the negative y-axis deriction.
                                  # the actuall programming process has nothing to do with this visualisation
                    A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                    return A_T
                elif p_a == d_c[1]: # refoemration parallelel to y-axis
                    A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                    return A_T
                elif p_a == d_c[2]: #leaves it be
                    return A
            elif rrank == 2:
                d_c = [[-1,0],
                       [0,-1]]
                if p_a == d_c[0]:
                    return A.T
                elif p_a == d_c[1]:
                    return A
            elif rrank == 1:
                return A

        e =  [1,2,3]
        if rank in e:
            return rotation_module(T, p_a)
        elif rank == 4:
            d_c  = [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]
            def rank4_rotation(T, t):
                def upward_shelve(T):
                    L = np.array([[T[i, j] for i in range(T.shape[0])] for j in range(T.shape[1])])
                    return L
                Tr = np.array([rotation_module(T[i], t) for i in range(T.shape[0])])
                L = upward_shelve(Tr)
                return L
            if p_a == d_c[0]:
                return rank4_rotation(T, [-1,0,0])
            elif p_a == d_c[1]:
                return rank4_rotation(T, [0,-1,0])
            elif p_a == d_c[2]:
                return rank4_rotation(T, [0,0,-1])
            elif p_a == d_c[3]:
                return T

    def Tensor_type(self, Tensor, args, margin=1.e-9, h=1.e-3):
        """
        (#ctv, #cv) = (contravariant, covariant)
        THis function take in a tensor and returns a tuple that says how
        many contravariant- and how many covariant components respectivelly
        the tesnor has.
        """
        X           = self.X
        J, Jt       = self.J, self.inv_J
        G           =  self.metric_tensor
        G_          =  self.inv_metric_tensor
        T           = Tensor
        rank        = np.ndim(T)
        d_X_ctv     = X(args + h) -X(args)
        d_B_ctv     = np.matmul(Jt, d_X_ctv)
        d_X_cv      = np.matmul(G_, d_X_ctv)
        d_B_cv      = np.matmul(J, d_X_cv)
        E = np.array([[d_X_ctv,d_X_cv,], [d_B_ctv, d_B_cv]])
        rank = np.ndim(T)

        def m(a,b, out=None):
            return np.matmul(a,b, out=None)
        def mtr(a,b, M=None, out=None):
            """
            M is an already formed matrix
            """
            if not M:
                return np.trace(m(a,b, out=None))
            elif not a and not b:
                return np.trace(M)

        if rank == 1:
            inv = m(T, d_B_cv) - m(m(J, T), d_X_cv)
            if np.abs(inv) < margin:
                return [0, 1]
            else:
                inv = m(T, d_B_ctv) - m(m(Jt, T), d_X_ctv)
                if np.abs(inv) < margin:
                    return [1,0]
                else:
                    raise Exception("""either the tensor in question is not a tensor or
                                    the math that the programmer implemented is flawed""")
        else:
            def co_contra_constructor(E, ind, rank):
                def help_module_rank2(a,ind):
                    """
                    ind equals either 1 or 0
                    """
                    if ind == 1 or ind ==0:
                        """
                        THisfunction creates a an array of matrices,
                        aka a 3D block of numbers that go into the mathematical teste function
                        """
                        _new_ = np.array([[a[ind,i]*a[ind,j].reshape(-1,1) for i in range(2)] for j in range(2)])
                        return _new_
                    else:
                        raise Exception("help_module_rank2: ind equals either 1 or 0")

                def help_module_rank3(b, ind):
                    """
                    This function creates an array of 3D blocks, t
                    his is a 4D array and boy will this be hard to work with, shit!!
                    """
                    if ind == 1 or ind ==0:
                        _new_ = help_module_rank2(b, ind)
                        __newer__ = np.array([[[_new_[i,j]*b[ind,k].reshape(-1,1) for k in range(_X_[i].shape[0])] for j in range(2)] for i in range(4)])
                        return __newer__
                    else:
                        raise Exception("help_module_rank3: ind equals either 1 or 0")

                if rank == 2:
                    return help_module_rank2(a, ind)
                elif rank==3:
                    return help_module_rank3(a, ind)

            def Mathematical_test_module(T1, T2, T3, rank):
                def scalar_project(T1, T2):
                     """
                     THis function projects rank 3 tensors onto each other and returns a scalar
                     This whole funcion makes me nervous
                     """
                     L1, L2 = T1.shape, T2.shape
                     m1 = np.array([[mtr(T1[i], T2[j]) for j in range(L2[0])] for i in range(L1[0])])
                     scalar_m1 = mtr(a=None, b=None, m1)
                     return scalar_m1

                def Tensor_permutator(A, rank):
                    a = A, A.T
                    if rank == 2:
                        L = np.array([[a[i],a[j]] for i in range(2) for j in range(2)])
                    elif rank == 3:
                        L = np.array([[a[i],a[j],a[k]] for i in range(2) for j in range(2) for k in range(2)])
                    elif rank == 4:
                        L = np.array([[a[i],a[j],a[k],a[p]] for i in range(2) for j in range(2) for k in range(2) for p in range(2)])
                    return L

                def TenMat_mult(T, M, R=None):
                    def single_tensor_single_matrix(T, M):
                        L_T = T.shape[0]
                        P =  np.array([m(T[i], M) for i in range(L_T)])
                        return P

                    def single_tensor_mult_matrix(T, R):
                        K = single_tensor_single_matrix
                        L_R = R.shape[0]
                        P = np.array([K(T, R[i]) for i in range(L_R)])

                    if not R:
                        return single_tensor_single_matrix(T, M)
                    elif not M:
                        return single_tensor_mult_matrix(T, R)

                if rank == 2:
                    R = Tensor_permutator(J, rank)
                    L = [mtr(T1, T2[i])- mtr(m(R[i,0], m(T, R[i,1])), T3[i]) for i in range(2**rank)]                       #VERY IMPORTANT
                    return L, R

                elif rank ==3:
                    R = Tensor_permutator(J, rank)
                    L = [scalar_project(T1,T2[i]) - scalar_project(TenMat_mult(T1, M=None, R=R[i]), T3[i]) for i in range(2**rank)]       #VERY IMPORTANT
                    return L, R

            def verify_tensor_type(r, matrix):
                """
                The way i'm using this fucntion makes it universal for any
                lenght of the tuple or list r.
                """
                E = [0,0]
                for i in r:
                    if (i-matrix < margin).all():
                        E[0] += 1
                    elif (i-matrix.T < margin).all():
                        E[1] += 1
                    else:
                        raise Exception("""something is awfully wrong with the code man!, this excetpion was raised in the rank 2 tensor
                                        type verification process where the tuple of the amuount of contra and covariant indicies was processed""")
                return E


            def Tensor_index_rotator(A, axis, lin_comb_axis=False, reverse=False):
                """
                LOL,This function is so good that you can just input a tensor twice and it will revert back
                to its original state, just like a matrix transpose; (M.T).T = M
                direction_list can only be [1,0,0], [0,1,0] or [0,0,-1] for rank=3,
                elsedirection_list =NOne
                """
                rank = np.ndim(A)
                direction_catalogue = [[-1,0,0],
                                       [0,-1,0],
                                       [0,0,-1]]
                K = L_A1,L_A2,L_A3 = A.shape[0], A.shape[1], A.shape[2]

                if rank == 2:
                    return A.T

                elif rank ==3:
                    if not lin_comb_axis: #this makes a false statement into a true statement
                        def transpose_machinery(A, direction_list):
                            if axis == direction_catalogue[-1]:
                                A_T = np.array([A[i].T for i in range(L_A1)])
                                return A_T
                            elif axis == direction_catalogue[-2]:
                                A_T = np.array([np.array([A[j].T[i] for j in range(L_A1)]).T for i in range(L_A3)])
                                return A_T
                            elif axis == direction_catalogue[-3]:
                                A_T = np.array([np.array([A[j,i] for j in range(L_A1)]).T for i in range(L_A2)])
                        return transpose_machinery(A, axis)
                    else:
                        if L_A1 == L_A2 == L_A3:
                            l = len(direction_catalogue)
                            L = np.array([transpose_machinery(A, direction_catalogue[k]) for k in range(l-1)])
                            LL = [np.zeros((L_A1,L_A1,L_A1)) + i for i in L]
                            return L[-1]
                        else:
                            raise Exception(f"your tensor must be a perfect cube,yours has the dimensions: {K}")


            def last_step(L, R, ind, rank):
                if L[ind] < margin:
                    tupple = verify_tensor_type(R[ind])                  #function call (4)
                    return tupple
                else:
                    raise Exception(f"""somehow the code didn't fuck-up anywhere and the main problem is therefore the math that
                                    i wanted to implement.This exception is in the last step of the rank {rank} tensortype verification
                                    process, the mathemaics is implemented in the sub-function with the name: Mathematical_test_module""")

            _B_ = co_contra_constructor(E, 1, rank)
            _X_ = co_contra_constructor(E, 0, rank)

            if rank == 2:
                T = T + T.T
                L, R = Mathematical_test_module(T, _B_, _X_)
                ind = L.index(min(L))
                return last_step(L, R, ind, rank)

            elif ranbk ==3:
                T = Tensor_index_rotator(T, None, lin_comb_axis=True)
                L, R = Mathematical_test_module(T, _B_, _X_)
                ind = L.index(min(L))
                return last_step(L, R, ind, rank)


    def TG(self,T, axis, L_R):
        """
        this function lowers or raises the index of a tensor of rank 1,2, or 3
        axis specifies which index and is a list
        L_R stands for Lower_Raise
        """
        rank = np.ndim(A)
        ttt = Tensor_type(self, T, args)
        p_a = axis
        def m(a,b, out=None):
            return np.matmul(a,b, out=None)

        TG_primer = self.TG_primer

        def Lower_Raise(T, p_a, L_R, ttt=ttt):
            def G_G(self, L_R):
                if L_R == [-1,1]:
                    return self.metric_tensor
                elif L_R == [1,-1]:
                    return self.inv_metric_tensor

            if L_R == [-1,1] or L_R == [1,-1]:  #----> lowering & raising indices respectivelly
                if ttt[1] == 0 and L_R[1] == -1 or ttt[0] == 0 and L_R[0] == -1:
                    raise Exception(f"you cant raise/lower an already raised/lowered tensor index. Your tensortype {ttt}")
                else:
                    T = TG_primer(T, p_a)   #------> this is the part where p_a is used

                    if rank == 1:
                        L = m(G_G(L_R), T)
                    elif rank == 2:
                        L = np.array([m(G_G(L_R), T)])
                    elif rank == 3:
                        L = np.array([m(G_G(L_R), T[i]) for i in range(T.shape[0])])

                    ttt = [ttt[0] +L_R[0], ttt[1] +L_R[1]]
                    return L,ttt
            else:
                raise Excperion("the following must hold: Lower_raise == [-1,1] or Lower_Riase == [1,-1]")
        return Lower_Raise(T, p_a, L_R)

    def direct_projection(self, T1, T2, p_a1, p_a2, contraction_order=1, p_a3=None):
        if contraction_order == 3 and not p_a3:
            raise Exception("you must specify p_a3 to be able to perform a self projection on a rank 4 tensor")
        else:
            def m(a,b, out=None):
                return np.matmul(a,b, out=None)
            def Gtr(a,b, M=None, out=None, M_direct=True):
                """
                M is an already formed matrix
                """
                if not M:
                    return np.trace(m(m(a,b), G))
                elif not a and not b:
                    return np.trace(m(M, G))
                elif direct:
                    return np.trace(M)

            r1 = np.ndim(T1)
            r2 = np.ndim(T2)
            def TG_primer(T, p_a, rank):
                    """
                    THis function primes the tensor for the lowering or the riase
                    of an index by the the function TG bellow. aka, it chooses
                    which index TG performs index juggling on.
                    """
                    def rotation_module(T, p_a, rank=rank):
                        A = T
                        if rank == 3:
                            d_c = [[-1,0,0],
                                   [0,-1,0],
                                   [0,0,-1]]
                            L_A0,L_A1,L_A2 = A.shape[0], A.shape[1], A.shape[2]
                            if p_a == d_c[0]: #the motivation for this is that: imagine ghat you had an axis in the negative x-direction attached to the cube
                                              # and thta you rotate clockwise, the face facing the forward direction will now
                                              #be facing the uppward direction, that is why the non-zero entries are -1.
                                              # the same logic applies for the next elif command but in the negative y-axis deriction.
                                              # the actuall programming process has nothing to do with this visualisation
                                A_T = np.array([np.array([A[j,i] for j in range(L_A0)]) for i in range(L_A1)])
                                return A_T
                            elif p_a == d_c[1]: # refoemration parallelel to y-axis
                                A_T = np.array([np.array([A[j].T[i] for j in range(L_A0)]) for i in range(L_A2)])
                                return A_T
                            elif p_a == d_c[2]: #leaves it be
                                return A
                        elif rank == 2:
                            d_c = [[-1,0],
                                   [0,-1]]
                            if p_a == d_c[0]:
                                return A.T
                            elif p_a == d_c[1]:
                                return A
                        elif rrank == 1:
                            return A

                    e =  [1,2,3]
                    if rank in e:
                        return rotation_module(T, p_a)
                    elif rank == 4:
                        d_c  = [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]
                        def rank4_rotation(T, t):
                            def upward_shelve(T):
                                L = np.array([[T[i, j] for i in range(T.shape[0])] for j in range(T.shape[1])])
                                return L
                            Tr = np.array([rotation_module(T[i], t) for i in range(T.shape[0])])
                            L = upward_shelve(Tr)
                            return L
                        if p_a == d_c[0]:
                            return rank4_rotation(T, [-1,0,0])
                        elif p_a == d_c[1]:
                            return rank4_rotation(T, [0,-1,0])
                        elif p_a == d_c[2]:
                            return rank4_rotation(T, [0,0,-1])
                        elif p_a == d_c[3]:
                            return T

            def TT_naive(T1, T2, p_a1, p_a2, r1, r2, contraction_order):
                T1 = TG_primer(T1, p_a1, r1)
                T2 = TG_primer(T2, p_a2, r2)

                e = [1,2]
                def MV_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order):
                    if contraction_order == 1:
                        T1 = TG_primer(T1, p_a, r1)
                        return m(T1, T2)
                    elif r1 != 2 and r2!= 2 and contrection_order == 2:
                        raise Exception("you cant contract twice if their first order contraction is a vector")
                    elif contrection_order == 2:
                        T1 = TG_primer(T1, p_a1)
                        T2 = TG_primer(T2, p_a2)
                        T3 =  m(T1, T2)
                        s  =  Gtr(a=None, b=None, M=T3)
                        return s

                def TT_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order, p_a3=p_a3):
                    def order4_trace(T1, T2, p_a3=p_a3):
                        L = np.array([[m(T1[i], T2[j]) for i in range(T1.shape[0])] for j in range(T2.shape[0])])
                        r3 = np.ndim(L)
                        L = TG_primer(L, p_a3, r3)
                        M = np.array([[Gtr(a=None, b=None, M=L[i,j]) for i in range(L.shape[0])] for j in range(L.shape[1])])
                        return M

                    if contraction_order==1:
                        if r2 == 1:
                            T1 = TG_primer(T1, p_a1, r1)
                            L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                            return L
                        else:
                            T2 = TG_primer(T2, p_a2, r2)
                            if r2 ==2:
                                L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                                return L
                            elif r2==3:
                                L = np.array([[m(T1[i], T2[j]) for i in range(T1.shape[0])] for j in range(T2.shape[0])])
                                return L
                    elif contraction_order == 2:
                        if r2 ==1:
                            T1 = TG_primer(T1, p_a1, r1)
                            L = np.array([m(T1[i], T2) for i in range(T1.shape[0])])
                            s = Gtr(a=None, b=None, M=L)
                            return s
                        elif r2 == 3:
                            return order4_trace(T1, T2)

                    elif contraction_order == 3:
                        a = order4_trace(T1, T2)
                        s = s = Gtr(a=None, b=None, M=a)


                if r1 in e and r1 in e:
                    return MV_comb(T1, T2, p_a1, p_a2, rank1, rank2)
                else:
                    return TT_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order)

            return TT_naive(T1, T2, p_a1, p_a2, r1, r2, contraction_order)

    def Christoffel_symbols(self, kind, args):
        G       = self.metric_tensor
        G_      = self.inv_metric_tensor
        D       = self.univ_finite_diff
        TG      = self.TG
        dirr    = self.direct_projection()

        Gamma_down  = 0.5*D(G, args)
        Gamma_up  = dirr(G_, Gamma1, [0,-1], [-1,0,0])

        if kind == 1:
            return Gamma_down
        elif kind == 2:
            return Gamma_up



    def Cov_deriv(self, T, args):
        pass


    def div(self, T, args):
        pass

    def curl(self, T, args):
        pass






























#-------------------------------------------------------------------------------------------------------------------------------------------

        #-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------------------






































































































"""

    def Tensor_contraction(self,Tensor1, Tensor2, index_reduction):
        G = self.metric_tensor
        G_ = self.inv_metric_tensor

        i_r = index_reduction
        a, b = Tensor1, Tensor2
        ra, rb = np.ndim(a), np.ndim(b)

        def matvec_comb(a,b, tensor_type_a, tensor_type_b):
            if
                b = np.matmul(G, b)
            elif tensor_type_b == "covriant" != tensor_type_a:
                b = npmatmul(G_,b)
            elif i_r == 1:
                if ra == rb ==1:
                    return np.matmul(a,b)
                if  ra == 2:
                    if rb == 1 or rb == 2:
                        return np.matmul(a,b)
            elif i_r == 2:

                if ra > 2 and rb > 2:
                    raise Exception("this function works with vaectors and matrices")
                elif ra ==2 != rb:
                    raise Exception("when index reduction==2, tensors should at least be matrices")
                else:
                    c = np.trace(np.matmul(a,b))
                    return c

        def rank3ten_matvec_comb(a,b):
            if i_r == 1:
                if tensor_type == "contravariant":
                    b = np.matmul(G, b)
                elif tensor_type == "covriant":
                    b = npmatmul(G_,b)

                l_a = a.shape[0]
                l_b = b.shape[0]
                if rb ==1 or rb ==2:
                    b = np.matmul(G, b)
                    L = np.array([np.matmul(a[i], b) for i in range(l_a)])
                elif rb == 3:
                    L = np.array([[np.matmul(a[i], np.matmul(G,b[j])) for i in range(l_a)] for j in range(l_b)])
                return L
            if i_r == 2:















        def rank4ten

        if i_r ==1:
            if ra == 2:
                return matvec_comb(a,b)
            if ra == 3:
                return tenmatvec_comb(a,b)

            if ra == 4:
                l_a = a.shape[0]
                if rb ==1 or rb == 2 or rb ==3:
                    L = np.array([])


        if i_r == 2:
            if ra < 2 and rb < 2:
                raise Exception("tensors need at least two indices to contract two indices")
            else:
                if ra ==2:
                    if rb == 2:
                        c = np.trace(np.matmul(a,b))
                        return c
                if ra ==3:
                    if r










    def Newt_motor(self, v0, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, finite=False, once=False):
            F = self.X
            J = self.J
            xvals = [v0]
            x_max = 1/margin
            k = 0

            def solving_motor(v):
                def finite_difference(v, h):
                    self.J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                                       [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                    return self.J

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








class function_Manifold(Manifold):
    def __init__(self, f, J, gradf, vec):
        self.vec  = vec
        self.f = f
        self.J = J
        self.gradf = gradf
        v = np.matmul(J, gradf)
        self.gradfp = v
        self.rmyu = np.vstack([J, v])

    def metric_tensor(self):
        Rm = self.rmyu

        G = np.matmul(Rm, Rm.T)
        return G
    def dim_f(self):
        pass



class Tensor(Manifold):
    def __init__(self, f, J, gradf, vec):
        super(). __init__(self, f, J, gradf, vec)





"""



















class Fractal2D_amer:
    """
    Task 1
    """
    def __init__(self, F, J):
        self.F = F
        self. J = J
        self.zero = []


    def Unified_PLOT(self, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT="both", A_method="int", coloroption="plasma"):
        """
        This Fucntion conatains all the other Tasks and methods  integrated into one function.
        """
        F = self.F
        J = self.J
        zero = self.zero

        def Newt_motor(self, v0, max_iter=max_iter, margin=1.e-4, step=step, h=1.e-9, finite=False, once=False):
            """
            Task 2
            """
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
                    self.J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                                       [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                    return self.J

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

        def Matrix_A(A_method):
            """
            Task 4
            """
            def tensor_motor(a,b,c,d,N):
                """
                Part of Task 4
                """
                x = np.linspace(a,b,N)
                y = np.linspace(c,d,N)
                X,Y = np.meshgrid(x,y)
                L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
                return  L

            if A_method == "int":
                def zero_collect_motor(F, J, v):
                    """
                    Task 2
                    """
                    v_a, v_i = Newt_motor(self, v)
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
                A = [[zero_collect_motor(L[i,j])[0] for i in range(N)] for j in range(N)]



            elif A_method == "abs":
                A = np.zeros((N,N))
                L = tensor_motor(a,b,c,d,N)
                for i in range(N):
                    for j in range(N):
                        new = Newt_motor(L[j,i])[0]
                        A[i,j] = np.abs(new[0]) + np.abs(new[1])
            return A

        if PLOT == True:
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, Matrix_A(A_method), cmap='plasma')

        elif PLOT == False:
            return Matrix_A(A_method)

        elif PLOT =="both":
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y,Matrix_A(A_method), cmap=f"{coloroption}")

            return zero




def Unified_PLOT(F, J, a,b,c,d,N, max_iter=100, margin=1.e-4, step=0.0001, h=1.e-9, approx=False, once=False, PLOT="both", A_method="int", coloroption="plasma"):
        """
        This Fucntion conatains all the other Tasks and methods  integrated into one function.
        """
        F = self.F
        J = self.J
        zero = self.zero

        def Newt_motor(F, J, v0, max_iter=max_iter, margin=1.e-4, step=step, h=1.e-9, finite=False, once=False):
            """
            Task 2
            """
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
                    self.J = np.array([[(F[0](v[0]+h,v[1]) -F[0](v))/h ,    (F[0](v[0],v[1]+h) -F[0](v))/h],
                                       [(F[1](v[0]+h,v[1]) -F[1](v))/h ,    (F[1](v[0],v[1]+h) -F[1](v))/h]])
                    return self.J

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

        def Matrix_A(A_method):
            """
            Task 4
            """
            def tensor_motor(a,b,c,d,N):
                """
                Part of Task 4
                """
                x = np.linspace(a,b,N)
                y = np.linspace(c,d,N)
                X,Y = np.meshgrid(x,y)
                L = np.array([[(X[i,j], Y[i,j]) for i in range(X.shape[0])] for j in range(X.shape[1])])
                return  L

            if A_method == "int":
                def zero_collect_motor(F, J, v):
                    """
                    Task 2
                    """
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
                A = [[zero_collect_motor(L[i,j])[0] for i in range(N)] for j in range(N)]



            elif A_method == "abs":
                A = np.zeros((N,N))
                L = tensor_motor(a,b,c,d,N)
                for i in range(N):
                    for j in range(N):
                        new = Newt_motor(L[j,i])[0]
                        A[i,j] = np.abs(new[0]) + np.abs(new[1])
            return A

        if PLOT == True:
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y, Matrix_A(A_method), cmap='plasma')

        elif PLOT == False:
            return Matrix_A(A_method)

        elif PLOT =="both":
            x = np.linspace(a,b,N)
            y = np.linspace(c,d,N)
            X,Y = np.meshgrid(x,y)
            fig, ax = plt.subplots()
            ax.pcolor(X, Y,Matrix_A(A_method), cmap=f"{coloroption}")

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
print(motor_integration(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both",A_method="int", coloroption="plasma"))
print(motor_integration(-2,2,-2,2,500, Test_F1, Test_J1, PLOT="both",A_method="int", coloroption="gray"))
print(motor_integration(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both",A_method="int", coloroption="plasma"))
print(motor_integration(-2,2,-2,2,500, Test_F2, Test_J2, PLOT="both",A_method="int", coloroption="gray"))
print(motor_integration(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both",A_method="int", coloroption="plasma"))
print(motor_integration(-2,2,-2,2,500, Test_F3, Test_J3, PLOT="both",A_method="int", coloroption="gray"))




print(motor_integration(-2,2,-2, 2,1000, Test_F, Test_J, PLOT="both",A_method="abs", coloroption="plasma"))
#print(motor_integration(-2,2,-2, 2, 10, Test_F3, Test_J3, PLOT="both",plotting_method="mine", coloroption="inferno_r"))
#print(motor_integration(-2,2,-2, 2,10, Test_F3, Test_J3, PLOT="both",plotting_method="fml", coloroption="gray"))






from scipy import *
import numpy as np
import matplotlib.pyplot as plt

plt.xlim=[-9000, 10000]
def f1(x,y):
    return np.array((x**3 - 3*x*y**2 -1, 3*x**2*y - y**3))

def j1(x,y):
    return np.array(((3*x**2-3*y**2, -6*x*y),
                   (6*x*y, 3*x**2-3*y**2)))

class fractal2D_erik:
    def __init__(self, F, J):
        self.f = f
        self.j = j
        self.zeroes = [] #our list of zeroes

    def newton(self, x0, maxiteration=1.e6, tol=1.e-4):
        x = np.array(x0, dtype=float) #makes it an array
        f_value=self.f(*x) #unzips the array and initialize it with our function
        f_norm = np.linalg.norm(f_value) #calculates our norm(length)
        iteration_counter = 0
        while abs(f_norm) > tol and iteration_counter < maxiteration:
            delta = np.linalg.solve(self.j(*x), -f_value) #solves for delta
            x += delta
            f_value = self.f(*x)
            f_norm = np.linalg.norm(f_value)
            iteration_counter += 1
        if abs(f_norm) > tol: #return -1 if the length is still not close to zero
            iteration_counter = -1
        return x #returns the coordinate

    def zeroes_index(self, x0, tol = 1.e-4):
        res = self.newton(x0)
        if (res==np.inf).any(): #check infinity
            return np.inf
        if not len(self.zeroes): #check if list is empty
           self.zeroes.append(res)
           return 0
        for k,i in enumerate(self.zeroes): #k is the iteration nr, i is the list value
            if np.linalg.norm(res-i) < tol: #check the lenght is less than tol
                return k
        self.zeroes.append(res)
        return len(self.zeroes) - 1 #needs to be -1 here because we append before the return

    def plot(self, N, a, b, c, d):
        x = np.linspace(a,b,N)#a,b defines the x-axis length and N how many points inbetween
        y = np.linspace(c,d,N)
        X, Y = np.meshgrid(x,y) #stores a grid in 2 matrices
        A = np.array([np.array([self.zeroes_index(np.array([X[i][j],Y[i][j]])) for i in range (N)]) for j in range(N)])
        #Stores a specific grid point which is our x0 of the form (X[i][j], Y[i][j])
        fig, ax = plt.subplots()
        ax.pcolor(X,Y,A)

    def simplified_newton(self, x0, maxiteration=1.e6, tol=1.e-4):
        pass
xa = np.array((1,5))
test1 = fractal2D(f1, j1)
test1.plot(100, -10,10,-10,10)
plt.show()


class OG_class(fractal2D_erik, Fractal2D_amer):
     def __init__(self,F, J):
          fractal2D_erik.__init__(self, f, j)
          Fractal2D_amer.__init__(self, F, J)


     def choose_method(self, a,b,c,d,N, amer=None, erik=None):
          if amer and not erik:
               return self.Unified_PLOT()
          elif erik and not amer:
               F = self.F
               self.plot(F, N, a,b,c,d)

















