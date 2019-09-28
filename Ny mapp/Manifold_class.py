
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


class Manifold:
    def __init__(self, X, J, args, H=None):
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
        self.args = args
        self.Hermitian = H
        self.X = X
        self.J = J
        Jt = J.T
        self.Jt = Jt
        G = np.matmul(J.T, J)(*X(*args))
        self.metric_tensor = G
        self.inv_metric_tensor = np.linalg.inv(G)
        self.detG = np.linalg.det(G)



    def Tensor_rank(self, T):
        rank = np.ndim(T)
        return rank

    def direction_catalogue(self,T):
        rank = self.Tensor_rank(T)
        d_c4  = [[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]

        d_c3  = [[-1,0,0],[0,-1,0],[0,0,-1]]

        d_c2 = [[-1,0],[0,-1]]
        if rank ==4:
            return d_c4
        elif rank ==3:
            return d_c3
        elif rank ==2:
            return d_c2


    def univ_finite_diff(self, Tensor, args, h=1.e-9):
        T = Tensor
        rank = np.ndim(T)

        def normal_deriv(Q):
            """
            This function calculates the usual partial derivative of a tensor valued function.
            It is divided into nested modules that, depending on the rank of the tensor-input,
            call deeper modules to do the approperiate job. The upper-most module turns vectors into jacobians of sorts,
            the second module uses the first module to turn the rows of a rank 2 tensor into
            the same type of jacobians and then turs them into a layered stack; a 3D-cube of derivates,
            and so on with higher order tensors.

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
                rows_m    = M.shape[0]

                Jacobi = np.array([[(A[i](*M[j]) - A[i](*args))/h  for j in range(rows_m)] for i in range(L_a)])
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
                return deriv_module1(Q)
            elif rank == 2:
                return deriv_module2(Q)
            elif rank == 3:
                return deriv_module3(Q)
            else:
                raise Exception(""""this class has no code for handling derivatives of
                                tensors of rank greater than 3, but it is easily fixable""")
        return normal_deriv(T)




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
                        delta = np.linalg.solve(univ_finite_diff(v, h), -F(v))
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

    def Tensor_type(self, Tensor, margin=1.e-9, h=1.e-3):
        """
        (#ctv, #cv) = (contravariant, covariant)
        THis function take in a tensor and returns a tuple that says how
        many contravariant- and how many covariant components respectivelly
        the tesnor has.
        """
        X           = self.X
        J, Jt       = self.J, self.inv_J
        G_          =  self.inv_metric_tensor
        T           = Tensor
        rank        = np.ndim(T)
        args        = np.array(args)
        d_X_ctv     = X(*(args + h)) -X(*args)
        d_B_ctv     = np.matmul(Jt, d_X_ctv)
        d_X_cv      = np.matmul(G_, d_X_ctv)
        d_B_cv      = np.matmul(J, d_X_cv)
        E = np.array([[d_X_ctv,d_X_cv,], [d_B_ctv, d_B_cv]])
        rank = np.ndim(T)

        def m(a,b, out=None):
            return np.matmul(a,b, out=None)
        def mtr(a,b, M=None, out=None):
            """
            M is (in its potential use) an already formed matrix
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

                    THisfunction creates a an array of matrices,
                    aka a 3D block of numbers that go into the mathematical test function down the road

                    """
                    if ind == 1 or ind ==0:

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
                    return help_module_rank2(E, ind)
                elif rank==3:
                    return help_module_rank3(E, ind)

            def Mathematical_test_module(T1, T2, T3, rank):
                def scalar_project(T1, T2):
                     """
                     THis function projects rank 3 tensors onto each other and returns a scalar
                     This whole funcion makes me nervous
                     """
                     L1, L2 = T1.shape, T2.shape
                     m1 = np.array([[mtr(T1[i], T2[j]) for j in range(L2[0])] for i in range(L1[0])])
                     scalar_m1 = mtr(a=None, b=None, M=m1)
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
                        return P

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
                            return LL[-1]
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

            elif rank ==3:
                T = Tensor_index_rotator(T, None, lin_comb_axis=True)
                L, R = Mathematical_test_module(T, _B_, _X_)
                ind = L.index(min(L))
                return last_step(L, R, ind, rank)


    def TG(self,T, axis, L_R, args_for_the_comparison_tensor_for_the_function_tensortype):
        """
        this function lowers or raises the index of a tensor of rank 1,2, or 3
        axis specifies which index and is a list
        L_R stands for Lower_Raise
        """
        args = args_for_the_comparison_tensor_for_the_function_tensortype
        rank = np.ndim(T)
        ttt = self.Tensor_type(T, args)
        p_a = axis
        G = self.metric_tensor(*args)
        def m(a,b, out=None):
            return np.matmul(a,b, out=None)

        TG_primer = self.TG_primer

        def Lower_Raise(T, p_a, L_R, ttt=ttt):
            def G_G(self, L_R):
                if L_R == [-1,1]:
                    return G
                elif L_R == [1,-1]:
                    return G

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
                raise Exception("the following must hold: Lower_raise == [-1,1] or Lower_Riase == [1,-1]")
        return Lower_Raise(T, p_a, L_R)

    def direct_projection(self, T1, T2, p_a1, p_a2, contraction_order=1, p_a3=None):
        G = self.metric_tensor(*args)
        if contraction_order == 3 and not p_a3:
            raise Exception("you must specify p_a3 to be able to perform a self projection (aka a trace) on a rank 4 tensor")
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
                elif M_direct:
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
                                              # the same logic applies for the next elif command but but now its in the negative y-axis' deriction.
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

                def MV_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order):
                    if contraction_order == 1:
                        T1 = TG_primer(T1, p_a1, r1)
                        return m(T1, T2)
                    elif contraction_order == 2:
                         if r1 != 2 and r2!= 2:
                             raise Exception("you cant contract twice if their first order contraction results in a scalar or a vector")
                         elif r1 == r2 == 2:
                             T1 = TG_primer(T1, p_a1)
                             T2 = TG_primer(T2, p_a2)
                             T3 =  m(T1, T2)
                             v  =  Gtr(a=None, b=None, M=T3)
                             return v

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

                e = [1,2]
                if r1 in e and r1 in e:
                    return MV_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order)
                else:
                    return TT_comb(T1, T2, p_a1, p_a2, r1, r2, contraction_order)

            return TT_naive(T1, T2, p_a1, p_a2, r1, r2, contraction_order)

    def Christoffel_symbols(self, kind, args):
        """
        kind = 1 gives you the fist kind
        kind = 2 gives you the second kind
        """
        J       = self.J
        G       = self.metric_tensor
        G_      = self.inv_metric_tensor
        D       = self.univ_finite_diff
        dirr    = self.direct_projection()
        H       = self.Hermitian
        d_c     = self.direction_catalogue


        if not H:
            Gamma_down   = 0.5*D(G, args)
            Gamma_up     = dirr(G_, Gamma_down, [0,-1], [-1,0,0])

            if kind == 1:
                return Gamma_down
            elif kind == 2:
                return Gamma_up
        else:
            Gamma_down   = dirr(H, J, d_c(H)[-1], d_c(J)[-1])
            Gamma_up     = dirr(G_, Gamma_down, [0,-1], [-1,0,0])

            if kind == 1:
                return Gamma_down
            elif kind ==2:
                return Gamma_up




    def Cov_deriv(self, T, args):
        """
        This function takes the  covariant derivative of a tensor of maximum rank 3

        """
        n,m         = self.Tensor_type(T, args)
        Gamma       = self.Christoffel_symbols(2, args)
        D           = self.univ_finite_diff
        dirr        = self.direct_projection
        rank        = np.ndim(T)
        d_c         = self.direction_catalogue
        dc_T        = d_c(T)
        dc_G        = d_c(Gamma)
        D_T         = D(T, args)

        if rank == 1:
            if n==1:
                T_Gamma = dirr(T, Gamma, None, dc_G[1])
            elif n==0:
                T_Gamma = -dirr(T, Gamma, None, dc_G[-1])
            return D_T + T_Gamma

        else:
            A = np.array([dirr(T, Gamma, dc_T[i], dc_G[2] for i in range(n-1))])
            B = np.array([dirr(T, Gamma, dc_T[i], dc_G[3] for i in range(n, m+1))])
            A.sum(axis=0)
            B.sum(axis=0)
            C = A+B
            return C


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



    def levi_civita(self, rank, dim, contravariant=None, covariant=None):
        """
        This fucntion creates the so called levi-civita tensor or epsilon permutation symbol.
        you can create it up to rank 4. The problem is only that i (the programmer) have mever seen
        how it looks beyond rank 3. for one, the rank 1 and rank 2 versions are good but according to wikipedia
        the rank 3 version should have only two numbers on every matrix in the cube.
        """

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











    def div(self, T, args):
        self.plopp = T
        pass

    def curl(self, T, args):
        pass

