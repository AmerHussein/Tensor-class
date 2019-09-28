"""
Created on Wed Jun 12 02:12:27 2019

@author: User
"""

from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import scipy.linalg as sl


def bisec(f, x1, x2, n, margin=1.e-11):
     if f(x1)<0 and f(x2)>0 or f(x1)>0 and f(x2)<0:
        for i in range(n):
            z = (x1 + x2)/2
            v = f(x1) * f(z)
            if v > 0:
                x1 = z
            elif v < 0:
                x2 = z
            delta = np.abs(x1 - x2)
            if delta < margin:
                return z
     else:
         raise Exception(f"f has the same sign at interval-endpoints, may not contain zero")





print(bisec(lambda x: x**3 +1, 2, 3, 10000))


























"""
def bisec(f,a,b,n, margin= 1.e-9):
    x=[a], y=[b]
    z = (a+b)/2
    if f(a)*f(b)>0:
        raise Exception("f may not be zero in interval [a,b]")
    elif f(a)*f(b) < 0:
        x.append(z)
        y.append(z)
        for i in range(n):
            F = f(x[-1])
            G = f(y[-1])
            if F < 0:
                x_new = (x[-1]+a)/2
                x.append(x_new)
            elif F > 0:
                a = x[-1]

            if G < 0:
                y_new = (y[-1] + b)/2
                y.append(y_new)
            elif G > 0:
                b = y[-1]

            if np.abs(f(x[-1])) < margin:
                if np.abs(f(y[-1])) < margin:
                    break
                    return [(x[-1], x[-2]), (y[-1], y[-2])]
                else:
                    return y[-1], y[-2]
            else:
                return x[-1], x[-2]n

print(bisec(np.cos, 0, np.pi/2,100))
"""



def power_decomposition(a):
    """
    This function takes in a number a and returns the part of the number
    that can be decomposed into a linear combination of powers of 10.

    ON PARAMETERS:
        a (float)

    ON RETURN:
        power_decomp (list):
            this list contains elements whose index values are the same as the
            exponent value that 10 is suposed to raised to.

    EXAMPLE:
        a = 1334.976
        resdiue(a)
        >>> [4,3,3,1]


    ATTENTION!:
        This function does not take into account the decimals of the input. (this was intentional yhough)

        THIS FUNCTION IS VERY BAD WHEN DEALING WITH LARGE NUMBERS AND NUMBERS
        THAT END WITH MORE THAN TWO ZEROS FOR SOME REASON , IT HAS SOME MAJOR FLAWS AND BUGS THAT
        MAKE THE ANSWER RESULT WEIRD (THE AMOUNT OF DIGITS SUDDENLY BECOME MORE THAN THE AMOUNT OF
        DIGITS OF THE INPUT AND SOME DIGITS MAY BE REPLACED BY "RADNOM" DIGITS) . BUT IT WORKS FOR SMALLER
        NUMBERS WITH LESS THAN 10 DIGITS.

    """
    K = np.array([10**(-j) for j in range(40)])
    L = np.array([10**(j) for j in range(40)])

    power_index  = np.where(a > L)
    power_value = len(list(power_index[0]))
    power_decomp =[]
    absa = np.abs(a)
    for i in range(1,power_value+1):
        for j in range(1000):
            absa -= 10**(power_value-i)
            if absa < 0:
                power_decomp.append(j)
                absa = absa + 10**(power_value-i)
                break
#            if absa == 0:
  #              power_decomp.append(0)
#        if  1<= np.abs(absa) <= 2:
 #           power_decomp.append(1)

    l = len(power_decomp)
    return power_decomp[::-1], l, power_value

A = 8990090987653463
print(power_decomposition(A))

def reverse_decomposition(A):
        """
        this function takes a list and uses its elements (numbers)
        and uses them as wheights to be multiplied with 10**(i) where i is the
        element's index inside the list.
        (it basically reverses the procedure of the function called: "power_decomposition",
        and returns a but with the decimals subtrracted away)

        ON INPUT:
            A (list)

        ON RETURN:
            val (int)
        """
        if type(A) != list:
            raise Exception("The object to be reverse decomposed should be a list")
        l = len(A)
        V = A
        val = 0
        for i in range(l):
            val += V[i]*10**(i)
        return val

B = reverse_decomposition([3, 6, 4, 3, 5, 6, 7, 8, 9, 0, 9, 0, 0, 9, 9, 8])

print(B)
print()

#def f(a,x):
#    return a*x
#x = np.linspace(0,1,1000)
#for i in range(10):
#    figure()
#    plot(x, f(i, x))
#    xlabel = ("x")
#    ylabel = ("f(i,x)")



from matplotlib.pyplot import *
from scipy import *
import scipy.integrate as integrate
import matplotlib.pyplot as plt
import numpy as np
from sympy import *
import scipy.linalg as sl
import numpy
from sympy.abc import x, y
from sympy import symbol, diff

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







def f(x):
    return -4*(1/np.sqrt(x**2 -4)) + x*(np.sqrt(x**2 -4) -1)/np.sqrt(x**2-4)



def fp(x, h=0.000001):
    deriv = (f(x+h) - f(x))/h
    return deriv
print(newton(f,fp , 5))


































def Newton(f,fp, x0, max_iter =400, div_deterrent =10**(10), margin =1.e-9):
    x =[x0]
    for i in range(max_iter):
        x_new = x[-1] -f(x[-1])/fp(x[-1])
        x.append(x_new)
        delta = np.abs(x[-1] - x[-2])
        travel_len = np.abs(x[0] - x[-1])
        if travel_len < div_deterrent:
            if i < max_iter:
                if delta < margin:
                    return x[-1], "convergent, x[-1] < margin"
            elif i == max_iter:
                if delta > margin:
                    return x[-1], "divergent, x[-1] > margin"
        elif travel_len > div_deterrent:

            raise Exception("x-values are travelling too far away")



























"""
LECTURE DAY 4

FUCNTIONS:
    - mathematical function:
        y = f(x):
            - x acan be a single variable or a multidimensional variable (such as a vecto)
            - x is called the argument/variable
            - y is called the dependent variable

    - the computer scinece way of defining a function:
    -    def parabola(x):
            return x**2

        - more pedeagogical definition
            def parabola(x):
                y = x**2
                return y
    - the input is also called "parameter"
    - to write parabola(0.3) is refered to as "calling the function"

    NOTE:
        in mathematics, we write "let f be a function" and never
        "let f(x) be a function". the former is the correct  and formal
        way.

    functions without arguments:
        def silly():
            return 3

        print(silly)   # 3 won't be returned
        print(silly()) # 3 will be returned

    Scope of  a varaible:
        - variables defined inside the definition of the function are said to belong to the
          scope of the function. rthey are unkown outside the function.
        -
          a = 3
          def f(x):
              x = 2
              return x

          if a variable is assigned before the definition of the function the value of that variable
          is knwon inside the fu ction defined afterwards. that variable is owned by the environment.
          of a function has an assignment of its parameter to a float, then when the environmental variable is
          is put through the function , value of the assignment isnide the function definition becomes
          the value of f(a). aka,
          a = 3
          a = x
          x = 2
          a = 2


    Default arguments:
        - if you have a default argument, other arguments should never have the same value as that
          default value.



    Docstring:
        - a d
        "sphinx" can be used to integrate latex text in your doctstrings












"""


help(plot)

def silly():
    return 3

z = silly()
print(z)

a = silly
w = a**2
print(a)
















"""
HOMEWORK 1 lOG
"""
def log_approx(x, n, diff = "on"):
    a,g = [(1+x)/2], [np.sqrt(x)]

    for i in range(n):
        a_new = (a[-1] + g[-1])/2
        a.append(a_new)
        g_new = np.sqrt(g[-1]*a[-1])
        g.append(g_new)
    if diff == "off":
        return (x-1)/ a[-1]
    elif diff == "on":
        return (x-1)/a[-1] -np.log(x)

x = np.linspace(2,5,100)



















































def sym_dif(a,b,swich="off"):
    if type(a) and type(b) != set:
        a,b = set(a), set(b)
    if swich =="off":
        sd_off = (a-b).union(b-a)
        return sd_off, "sd_off"
    elif swich == "on":
        sd_on = a.symmetric_difference(b)
        return sd_on, "sd_on"

a = set([1,2,3,4])
b = set([3,4,5,6])
print(sym_dif(a,b))
print(a.intersection_update(b))
print(a.intersection(b))



"""
def dup_red(s, d=[]):
    diag = [0]
    if s not in d and s not in diag:
        d.append(s)
    return np.array(d)

#vdup_red= np.vectorize(dup_red)

A = np.array([0,20,30,40])
print(dup_red(3))

vdup_red = np.vectorize(dup_red)
print(list(vdup_red(A)))
"""






































def Beta_bool(B,n, max_iter =30, margin = 1.e-9, x0=1):
    """
    denna funktion genomför n rekusrion under en for-loop
    ON INPUT:
        B (float)
    """
    c= [0.2  ,5]
    x=[x0]
    for i in range(n):
        x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
        x.append(x_new)
        delta = np.abs(x[-1]- x[-2])
        if i < max_iter:
            if delta < margin:
                return True, x

        elif i == max_iter:
            if delta > margin:
                return False, x

print(Beta_bool(0.05, 100))


def Beta(B,n, p = "x", margin = 1.e-9, x0=1):
    """
    this function executes a recursion with the aid of a for-loop
    ON INPUT:
        *B (float):
            it is very important that we keep the star.
    ON OUTPUT:
        x (list)
        or
        x[-1] (float)
        or
        conv: {x[-1]} (string)

    """
    c,x = [0.2  ,5], [x0]
    for i in range(n):
        x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
        x.append(x_new)
        delta = np.abs(x[-1]- x[-2])
        if delta < margin:
            break

    if p == "x_lastval_sring":
        return f"""conv: {x[-1]}"""

    elif p == "x_lastval":
        return x[-1]

    elif p == "x":
        return x
print(Beta(0.5, 100))



def sign_separation(A, neg=[], zero=[], pos=[]):
    for s in A:
        if s > 0:
            pos.append(s)
        elif s < 0:
            neg.append(s)
        else:
            zero.append(s)

    if len(zero) == 0:
        return [neg, pos]
    else:
        return [neg, zero, pos]

b = [-0.5, 0.5, -0.25, 0, 0.25]
print(sign_separation(b))





def seq_conv(n, *var, p = "x_lastval", sign_separation_swich ="off", bool_convert = "no", margin = 1.e-9, x0 =1):
    """
    Denna funtion tar en lista *var, och genomför funcktionen Beta över dess element
    ON INPUT:
        n (int):
            amount of itertions
        var (list)
        p (string)

    ON OUTPUT:
        bb (list)
    """

    bb = []
    for s in var:
        s = Beta(s,n)
        bb.append(s)

    if bool_convert == "yes":
        V = []
        for s in var:
            v = Beta_bool(s,n)
            V.append(v)
        return


    if sign_separation_swich == "on":
            return sign_separation(bb)
    elif sign_separation_swich == "off":
        return bb

b = [-0.5, 0.5, -0.25, 0, 0.25]
print(seq_conv(100, *b))








"""
LECTURE NOTES day 3
last lecture notes were deleted by an accidental shut down.(forgot to charge)
CONTAINER TYPES:
    SLICING:
        - syntax:
            - L = [’C’, ’l’, ’o’, ’u’, ’d’]
            type(L) #list
            L[i:j]

    DICTIONARIES:
        - a dictionary is an undordered structure if ( key, value)
          pair. ONe addresses elemts in a dictionary with their respective key
        - the item can be any object
        - the keys can only be immutable objects, such as strings or  tuples
          example:
              homework = {"Anna":  True, "kersin": False}

              homework["Anna"]    # returns the value/item associated with th key "Anna": True
              homewrk["kerstin"]  # returns False
              - You cana also delete value/item by writing:
                  del homework["Anna"]
        - to say "what item comes first makes no sense, because the items are orderd not by
          index but rather by keys.

        TYPICAL dictionary METHODS: items, keys, values:
            - d={’name ’:’Elsa ’,’age ’:23 ,’hobby ’:’fishing ’}
              "   print(f"a typical student: {list(d.values())})   "
              - by suing a for-loop:
                  for value in d.values():
                      print(f"a typical student: {values}")

         - dictionaries are mainly used for providing functions with arguments in a compact way
           (see Unit 4).
         - to collect options of a method: options={’tol’:1.e-3,’stepsize’:0.1,’maxit’:1000}

    SETS:
        A set is a colection f well defined and distinct objects, considered as an object in its own roght.

        - fruitbasket = set(["apple", "pear", ])
        - The most important operation is in, meaning ∈ (is element of):
            ’plum ’ in fruitbasket # returns False
            ’pear ’ in fruitbasket # returns True

    TUPLES: you cannot change elements in a











STRING FORMATTING:
    a = "{a}\n b=3" is the same as two strings
"""

"""








"""































































"""
def sign_separation(s,neg=[], zero=[], pos=[]):

    if s > 0:
        pos.append(s)
    elif s < 0:
        neg.append(s)
    else:
        zero.append(s)

    neg = np.array(neg)
    zero = np.array(zero)
    pos = np.array(pos)

    return np.array([[neg], [zero], [pos]])
vsign_separation = np.vectorize(sign_separation)


B = np.array([-0.5, 0.5, -0.25, 0.25])
#print(vsign_separation(B))

A = vsign_separation(B)
C = list(A)

print("---------------------")
print(C[0])
print("---------------------")

print(C[2])
"""









































def seq_conv(n, *var, p = "x_lastval", sign_separation_swich ="off", iteration_regulation = "on", margin = 1.e-9, x0 =1):
    """
    Denna funtion tar en lista *var, och genomför funcktionen Beta över dess element

    ON INPUT:
        n (int):
            amount of itertions
        var (list)
        p (string)

    ON OUTPUT:
        bb (list)

    """
    x,c,b = [x0], [0.2  ,5], []




    if iteration_regulation == "on":
        def Beta_moderator(B, max_iter = 30):
            """
            denna funktion genomför n rekusrion under en for-loop
            ON INPUT:
                B (float)
            """
            for i in range(max_iter+1):
                x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
                x.append(x_new)
                if np.abs(x[-1]- x[-2])< margin:
                    break
                    return True
                else:
                    return False

    elif   iteration_regulation == "off":
        def Beta(B):
            """
            denna funktion genomför n rekusrion under en for-loop
            ON INPUT:
                B (float)
            """
            for i in range(n):
                x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
                x.append(x_new)
                if np.abs(x[-1]- x[-2])< margin:
                    break

            if p == "x_lastval_sring":
                return f"""conv: {x[-1]}"""

            elif p == "x_lastval":
                return x[-1]

            elif p == "x":
                return x

    """
    here we apply Beta to our input var:
    """
    for s in var:
        s = Beta(s)
        b.append(s)




    if sign_separation_swich == "yes":
        """
        this functionseparates the elemts of list, A, into a list containing three lists,
        with negative, zero elements and positive elemts respectivelly
        """
        def sign_separation(A, neg=[], zero=[], pos=[]):
            for s in A:
                if s > 0:
                    pos.append(s)
                elif s < 0:
                    neg.append(s)
                else:
                    zero.append(s)

                if len(zero) == 0:
                    return [neg, pos]
                else:
                    return [neg, zero, pos]
        return sign_separation(b)

    elif sign_separation_swich == "no":
        return b

B = [-0.5, 0.5, -0.25, 0.25]
print(seq_conv(100, *B))

"""
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
__________________________________________________________________________________________________________________________________________
"""

"""
task  4, 5,7
"""

def sign_separation(*A, neg=[], zero=[], pos=[]):

    for s in A:
        if s > 0:
            pos.append(s)
        elif s < 0:
            neg.append(s)
        else:
            zero.append(s)

    if len(zero) == 0:
        return [neg, pos]
    else:
        return [neg, zero, pos]



def Beta(B,n, margin, x0):

    """
    this function executes a recursion with the aid of a for-loop
    ON INPUT:
        *B (float):
            it is very important that we keep the star.
    ON OUTPUT:
        x (list)
        or
        x[-1] (float)
        or
        conv: {x[-1]} (string)

    """
    c,x = [0.2  ,5], [x0]
    for i in range(n):
        x_new = c[0] * x[-1] -B*(x[-1]**2 - c[1])
        x.append(x_new)
        if np.abs(x[-1]- x[-2])< margin:
            break

    if p == "x_lastval_sring":
        return f"""conv: {x[-1]}"""

    elif p == "x_lastval":
        return x[-1]

    elif p == "x":
        return x



def seq_conv(n, *var, p = "x_lastval", sign_separation_swich ="on", iteration_moderator = "off", margin = 1.e-9, x0 =1):
    """
    Denna funtion tar en lista *var, och genomför funcktionen Beta över dess element
    ON INPUT:
        n (int):
            amount of itertions
        var (list)
        p (string)

    ON OUTPUT:
        bb (list)
    """
    bb = []
    for s in var:
        s = Beta(s)
        bb.append(s)

    if sign_separation_swich == "on":
            return sign_separation(bb)
    elif sign_separation_swich == "off":
        return bb






"""
task 3
"""


def iter_(f, max_iter, margin =1.e-9):
    x = [f(0)]
    for i in range(1,max_iter):
        x_new = f(i)
        x.append(x_new)
        if np.abs(x[-1]) < margin:
            break
    return x, min(x), len(x)


def g(x):
    return ((np.sin(x))**2)/x
print(iter_(g, 1000))





"""
exersise 1
"""


"""
task 1
"""


def fix_point(f, x0, max_iter,  margin = 1.e-10, MARGIN = 1.e-7, min_iter =30):


    """
    This function takes a function and an initiat guess on the the first x-value and returns a fixed point
    with the property  x = f(x)
    """


    N = max_iter*2
    x = [x0]
    if c == "index":
        for i in range(max_iter):
            x.append((f(x[i-1])))
            delta = np.abs(x[-1] -x[-2])
            if delta <= margin and i > min_iter:
                break
            return x[-1], delta, f"""delta < {margin}""".format(margin=margin)

    elif c=="claus_method":
        for i in range(max_iter):
            x_new = f(x[-1])
            x.append(x_new)


    DELTA = np.abs(x[-1] -x[-2])
    if  margin <= DELTA <= MARGIN:
        for i in range(max_iter, N):
            x.append((f(x[i-1])))
            D = np.abs(x[-1] -x[-2])
            if D < margin and i > min_iter:
                break
            return x[-1], """fixed DELTA to delta""" , D

    elif DELTA >= MARGIN:
        return x[-1], DELTA, "DELTA"
    else:
        D = np.abs(x[-1] -x[-2])
        return x, D, "D"


def f(x, a = 0.5):
    return np.sin(x) - a*x +30

print(fix_point(f, 0.5, 300))























"""
sommarkurs i python: NOTES [computational programming]:
    - [(8:15 - 10   Lecture), (10.15 - 12teaching Assistant)]  = unit
    - OBS! på onsdag så börjar lektionen kl 8 och slutar 9
    - HOMEWORK (mandatory) (group of 2) - present to the teaching assistance
    - FINAL PROJECT (madatory) - organise a code in a group (2-3)
                             - presen to the proffessor

    def ?(question):
        return "claus answer next lecture"

    - programming will be stongly related to mathematcial computing


Spyder:
    - Spyder is what we call a work bench
    - These notes that (you) read now reside in the "editor"
    - the box ---> is called the "sandbox" (console):
        - the condole acts as the command promt. you can write lists and
          write mathematical operations

    - command window. black window, in windows power shell,  write commands directly:
        -write python commands
        - Write "powershell" in the search bar down left
        - and inside the command promt (the balck window) write: ipython

    - red cikles: code will never work
    - yelow triangles: gives you hints to make code better

Objects:
    - numbers:
        - floating points (reals)
        - complex flotating point (complex)
        - integers
    - Arrays:
        - vectors
        - matrices
    - Lists:
        - can contain objects of different types

        - List utilities:
            - the function "range(a, b)" can generate a list of numbers from a to b
            - List comprehension:
                example:
                    L = [1,2,3,4]
                    L2 = [x*2 for x in L]  # [2,4,6,8]
                    L3 = [x for in L if 2 < x <= 10]  #[3,4]
            - Concatenation:
                A = [1,2]
                B = [4,5]
                L = A + B    # [1,2,4,5]
                K = A * 3    # [1,2,1,2,1,2]
            - For loop:
                L = [1,2,3,4]
                for s in L:
                    print(s*2)    # prints: 2,4,6,8

                n = 30
                for i in range(30):
                    do_something
Indentation (amer you already know)

Basic Plotting:
    - We generate two lists:
        x_list = list(range(100)) #the first 100 number
        y_list = [sqrt(x) for x in x_list]

    - And then we make a grpah:
        plot(x_list, y_list)
        title("my first plot")
        xlabel("x")
        ylabel("sqrt(x))
        show()

"""




def TASK(number, *args):

    """
    Task 1
    """
    def is_zero(f, a, tol = 1.e-10):
        if f(a) == 0:
            return True
        elif f(a) != 0 and np.abs(f(a)) < tol:
            return "value is very close to zero"
        else:
            return False

    """
    Task 6,8,9
    """
    def equidist_plot_y(f, a,b, n, c="linspace", Y = "array"):
        """
        This function creates two lists, or two arrays depedning on the value of
        Y.
        ON INPUTS:
            f (function):
                mathematical function
            a (float):
                lower interval bound
            b (float):
                upper interval bound
            n (float):
                array or "list" dimension
            c (string):
                can have the following values:
                    linspace, partition, counter
                is set by default as linspace
            Y (string):
                can have the following values:
                    array, list
                is set by default as array

        ON RETURN:
            depends on Y
        """

        if c=="linspace":
            L = list(np.linspace(a,b, n))
            u = np.array(L)

        elif c=="partition":
            L = [a+(b/n)*k for k in range(n+1)]
            u = np.array(L)

        elif c== "counter":
            counter, L = a, []
            for k in range(n):
                counter += (b/n)*k
                L.append(counter)
                u = np.array(L)

        vf = np.vectorize(f)
        y = vf(u)
        if Y == "array":
            u = list(u)
            y = list(y)
            return np.array([u,y])
        elif Y == "list":
            return [list(u), list(y)]

    def sum(m,n, f):
        if n > m:
            init = f(m)
            Q = np.zeros(1+n-m)
            for i in range(m,n):
                init += f(i)
                Q[i] = init
            return Q[-1], Q
        else:
            raise Exception("n should be greater than m")

    """
    Task 12
    """
    def list_gen(n, h =1000, a = -0.5, c ="append", d = None):
        """
        This function returns an array or a list depending on the value of d
        """
        u = [1, np.exp(h*a), np.exp(2*h*a)]
        td = [x*h for x in range(n)]
        Td = np.exp(a*np.array(td))


        """
        three ways creating u
        """
        if c == "append":
            for i in range(3,n):
                u.append(u[i-1] + h*a*( (23/12)*u[i-1] - (4/3)*u[i-2] + (5/12)*u[i-3] ))

        elif c == "index":
            u = np.zeros(n)
            for i in range(3,n):
                u[i] = u[i-1] + h*a*( (23/12)*u[i-1] - (4/3)*u[i-2] + (5/12)*u[i-3] )

        elif c == "claus_method":
            B = [23/12,4/3, 5/12]
            for i in range(n):
                u_new = u[-1] + h*a*(B[0]*u[-1] - B[1]*u[-2] + B[2] *u[-3])
                u.append(u_new)


        """
        creating td and Td
        """
        M = np.array([[Td], [u]])
        delta = list(np.abs(M[1] -M[0]))


        """
        the dependence of return on the value of d is formulated using if/elif statemts
        """
        if d == "u":
            return u
        elif d == "Td":
            return Td
        elif d =="td":
            return td
        elif d == "M":
            return M
        elif d == "delta":
            return delta
        elif d == "pure confusion":
            return np.array([u, Td, td, M, delta])

        """
        returning results depending on what number the task is.
        """

    if number == "1":
        return is_zero(*args)

    elif number =="6" or "8" or "9":
        return equidist_plot_y(*args)
    elif number == "12":
        return list_gen(*args)





















def fix_point(f, x0, max_iter,  margin = 1.e-10, MARGIN = 1.e-7, min_iter =30):

    """
    This function takes a function and an initiat guess on the the first x-value and returns a fixed point
    with the property  x = f(x)
    """
    N = max_iter*2
    x = [x0]
    if c == "index":
        for i in range(max_iter):
            x.append((f(x[i-1])))
            delta = np.abs(x[-1] -x[-2])
            if delta <= margin and i > min_iter:
                break
            return x[-1], delta, f"""delta < {margin}""".format(margin=margin)
    elif c=="claus_method":
        for i in range(max_iter):
            x_new = f(x[-1])
            x.append(x_new)


    DELTA = np.abs(x[-1] -x[-2])
    if  margin <= DELTA <= MARGIN:
        for i in range(max_iter, N):
            x.append((f(x[i-1])))
            D = np.abs(x[-1] -x[-2])
            if D < margin and i > min_iter:
                break
            return x[-1], """fixed DELTA to delta""" , D

    elif DELTA >= MARGIN:
        return x[-1], DELTA, "DELTA"
    else:
        D = np.abs(x[-1] -x[-2])
        return x, D, "D"

    """
    else:
        l = len(x)
        D = np.zeros(l)
        for i in range(1, l-5):
            D[i] = np.abs(x[-i] -x[-(i+1)])
            if D[i] < margin:
                break
                DD = list(D[1:i+1])
                min_error = min(DD)
                return x, min_error, "min_error"
    """




def f(x, a = 0.5):
    return np.sin(x) - a*x +30

print(fix_point(f, 0.5, 300))















"""
-------------------------------------------------------------
"""


"""
Task 12
"""
def list_gen(n, h =1000, a = -0.5, c ="append", d = None):
    """
    This function returns an array or a list depending on the value of d
    """
    u = [1, np.exp(h*a), np.exp(2*h*a)]
    td = [x*h for x in range(n)]
    Td = np.exp(a*np.array(td))


    """
    two ways creating u
    """
    if c == "append":
        for i in range(3,n):
            u.append(u[i-1] + h*a*( (23/12)*u[i-1] - (4/3)*u[i-2] + (5/12)*u[i-3] ))

    elif c == "index":
        u = np.zeros(n)
        for i in range(3,n):
            u[i] = u[i-1] + h*a*( (23/12)*u[i-1] - (4/3)*u[i-2] + (5/12)*u[i-3] )

    """
    creating td and Td
    """
    M = np.array([[Td], [u]])
    delta = list(np.abs(M[1] -M[0]))


    """
    the dependence of return on the value of d is formulated using if/elif statemts
    """
    if d == "u":
        return u
    elif d == "Td":
        return Td
    elif d =="td":
        return td
    elif d == "M":
        return M
    elif d == "delta":
        return delta
    elif d == "pure confusion":
        return np.array([u, Td, td, M, delta])

print(list_gen(5, h =1000, a = -0.5, c ="append", d="pure confusion"))

"""
----------------------------------------------------------------
"""


"""
task 1
"""
def is_zero(f, a, tol = 1.e-10):
    if f(a) == 0:
        return True
    elif f(a) != 0 and np.abs(f(a)) < tol:
        return "value is very close to zero"
    else:
        return False

f = lambda x: x**2 + 0.25*x -5
print(is_zero(f, 2.3))



"""
Task 2
"""

"""
----------------------------------------------------------------
"""


"""
Task 6,8,9
"""
def equidist_plot_y(f, a,b, n, c="linspace", Y = "array"):

    """
    This function creates two lists, or two arrays depedning on the value of
    Y.
    ON INPUTS:
        f (function):
            mathematical function
        a (float):
            lower interval bound
        b (float):
            upper interval bound
        n (float):
            array or "list" dimension
        c (string):
            can have the following values:
                linspace, partition, counter
            is set by default as linspace
        Y (string):
            can have the following values:
                array, list
            is set by default as array

        ON RETURN:
            depends on Y
    """

    if c=="linspace":
        L = list(np.linspace(a,b, n))
        u = np.array(L)

    elif c=="partition":
        L = [a+(b/n)*k for k in range(n+1)]
        u = np.array(L)

    elif c== "counter":
        counter, L = a, []
        for k in range(n):
            counter += (b/n)*k
            L.append(counter)
        u = np.array(L)


    vf = np.vectorize(f)
    y = vf(u)
    if Y == "array":
        u = list(u)
        y = list(y)
        return np.array([u,y])
    elif Y == "list":
        return [list(u), list(y)]

[x_list, y_list] = equidist_plot_y(np.arctan, 0,1, 100, c="partition")
print(type(x_list))
x = list(x_list)
print(type(x))

"""
----------------------------------------------------------------
"""




print(equidist_plot_y(np.arctan, 0,1, 100, c="partition", Y ="list"))


"""
Task 9
"""
"""
----------------------------------------------------------------
"""
def sum(m,n, f):
    if n > m:
        init = f(m)
        Q = np.zeros(1+n-m)
        for i in range(m,n):
            init += f(i)
            Q[i] = init
        return Q[-1], Q
    else:
        raise Exception("n should be greater than m")

f = lambda x: 1/np.sqrt(x)
print(sum(1,100, f))
"""
----------------------------------------------------------------
"""






"""
----------------------------------------------------------------
"""



def f(x):
    """
    vectorizing an already vectorized function: arctan
    """
    return np.arctan(x)
vf = np.vectorize(f)
a = np.array([1,2,3])
print(vf(a))















"""
----------------------------------------------------------------
"""

def zero(f,a,b,c="efficient",div= 1000, error_margin = 10**(-2)):
    """
    ----------------------
    denna funktion beräknar nollställen av en funtion som f.
    ON INPUTS:
        f (function):
            any mathematical function is not too crazy, of one variable.
        a (float):
            lower interval bound
        b (float):
            upper interval bound
    ON RETURN:
        let's call the wanted result "val", and the error "e"
        (val-e, val, val+e) (tuple):
            gives you the error bounds of the result.
    OBS!!:
        du kan välja mellan en snabb version och en långsam version. den snabba
        versionen använder sig av vectorization och list comprehension och annat kul i guess.
        c = "efficient":  så blir koden snabb
        c = "inefficient": så blir koden mcycket långsammare förstora värden på div.


    ---------------------
    """

    if np.abs(f(a)) <= error_margin:
        return a

    elif np.abs(f(b)) <= error_margin:
        return b

    elif c == "inefficient":
        x = list(np.linspace(a,b, div))
        func_värden = list(np.zeros_like(x))
        r = []
        rr = []
        l = len(x)
        for k in range(l):
            func_värden[k]= f(x[k])
        for i in range(l):
            if np.abs(func_värden[i]) < error_margin:
                r.append(np.abs(func_värden[i]))
                rr.append(x[i])
        _min_ = min(r)
        ind = r.index(_min_)
        return (rr[ind]- error_margin, rr[ind] , rr[ind] + error_margin),"---", r

    elif c == "efficient":
        vf = np.vectorize(f)
        x = np.linspace(a,b,div)
        f_vals = vf(x)

        def find_lowest_value(x):
            if x < error_margin:
                return np.abs(x)
            else:
                return np.inf

        vfind_lowest_value = np.vectorize(find_lowest_value)
        Q = vfind_lowest_value(f_vals)
        q = list(Q)
        _min_ = min(q)
        ind = q.index(_min_)
        return (x[ind-1]  , x[ind], x[ind+1])



"""
----------------------------------------------------
"""










































"""
---------------------------------------------------------
"""






def linalg_op(A, B, C = "plagiat"):

    """
    --------------------------
    YOYO wassup my niqqa! SIT THE FUCK DOWN AND READ!!

    Linear algebra operations: dot prod, matrix-vector mult, matrix-matrix mult
    and rank 3 tensor-matrix mult. This is a binary operation...
    (i should work on extending it to multiple array-inputs)

    on inputs:
        A and B (np.array()) given arrays compatible with above function description.
        C (string):
            If you input the string "plagiat" the function, linalg_op,  will use the built in python NumPy
            linear algebra operations.

            If input the string "ejplagiat" the function, linalg_op, will use analogous operations
            built by Amer.

            OBS!:
                C is by efault set to be "plagiat" so that the operations are computed faster
                and may be changed by the user bys just input C = "ejplagiat" as a positional argument
    on return:
        Q (np.array() or float) the result is either an array or scalar depending on the input
    --------------------------
    """

    def vv(A, B):
        if type(A) != type(B) != np.array():
            A = np.array(A)
            B = np.array(B)
        elif len(A) != len(B):
                raise Exception("your input vectors should have the same dimesion")
        return sum(list(A * B))

    def mv(M,A):
        L = M * A
        l = len(A)
        e = np.ones(l)
        q = np.zeros(l)
        for i in range(l):
            q[i] = vv(L[i], e)
        return q

    def mm(M, W):
        L,P = M.shape[1], W.shape[0]
        Q = np.zeros((L,P))
        for i in range(P):
            Q[i] = mv(M,W[i])
        return Q

    def tm(T, M):
        L,P,P_ = T.shape[2], M.shape[0], M.shape[1]
        Q = np.zeros((L,P,P_))
        for i in range(L):
            Q[i] = mm(T[i], M)
        return Q

    if C == "plagiat":
        return np.matmul(A,B), help(linalg_op)

    elif C == "ejplagiat":
        if np.ndim(A) == 1 and np.ndim(B) ==1:
            return vv(A,B), help(linalg_op)

        elif np.ndim(A) == 1 and np.ndim(B) ==2:
            return mv(B,A), help(linalg_op)

        elif np.ndim(A) == 2 and np.ndim(B) ==1:
            return mv(A, B), help(linalg_op)

        elif np.ndim(A) == 2 and np.ndim(B) ==2:
            return mm(A,B), help(linalg_op)

        elif np.ndim(A) == 3 and np.ndim(B) ==2:
            return tm(A,B) , help(linalg_op)

"""
----------------------------------------------------------------------
"""







S = np.array([[0,3,5,11],
              [3,0,13,17],
              [5,13,0,27],
              [11,17,27,0]])
a = np.array([3,0,13,17])
print(linalg_op(S, a, C = "ejplagiat"))

I = np.identity(4)
a = np.array([3,0,13,17])
b = np.array([0,3,5,11])
print(vv(a,b))
print(I *a)
print(mv(S,a))










#x, y, z = symbols('x y z', real=True)
#f = 4*x*y + x*sin(z) + x**3 + z**8*y
#delxf = f.diff(z)
#print(delxf)


#*args = symbols('*args', real=True)
#f = 4*x*y + x*sin(z) + x**3 + z**8*y
#delxf = f.diff(z)
#print(delxf)












def add(*args):
    q =0
    for x in args:
        q += x
    return q
x = np.arange(4)
print(add(*x))






#---------------------------------------------------------------------
def rel(F):
    K = 2*shape(F)[0]
    a = int((K-4)/2+1)
    b = int(-(K-4)/2)
    e = np.identity(F.shape[0])*1.e-10
    B = (F-e >0)
    Q = np.array(F[B])



    return np.hstack([Q[:a],  Q[b:]]).reshape()

S = np.array([[0,3,5,11],
              [3,0,13,17],
              [5,13,0,27],
              [11,17,27,0]])

print(rel(S))
#----------------------------------------------------------------------


def hej(x,y):
    if type(x) != np.ndarray:
        if type(x) == list:
            x = np.array(x)
        else:
            raise Exception("first imput should be a vector")

    if type(y) != np.ndarray:
        if type(y) == list:
            y = np.array(y)
        else:
            raise Exception("second imput should be a vector")

    if x.size > y.size:
        i = np.ones(x.size - y.size)
        v = np.hstack([y,i])
        W = np.ones((len(x), len(v)))
        L = W*x +v.reshape(-1,1)
        LL = L[0:len(y),:]
        return LL

    elif y.size > x.size:
        i = np.ones(y.size - x.size)
        u = np.hstack([x,i])
        W = np.ones((len(u), len(y)))
        L = W*u +y.reshape(-1,1)
        LL = L[:, 0:len(x)]
        return LL

    elif y.size == x.size:
        W = np.ones((len(x), len(y)))
        L = W*x +y
        return L

r = np.array([0,1])
p = np.array([0,1,2])
print(hej(p,r) == hej(r,p))
print(hej(p,r))
print(hej(r,p))

    #AAAND of course , there has to be a way to fucking do this hit much faster and more reliably. fuck me

def hhej(x,y):
    return x.reshape(-1,1) + y

x = np.array([1,2,3,4])
y = np.array([2,3,4])
print(hej(y,x))



P = np.array([[0,3,5],
              [3,0,13],
              [5,-13,0],
              [-11,17,27]])
v = np.array([1,2,3,4])

G = P+v.reshape(-1,1)
print(v.reshape(-1,1))
print("---")
print(G)




def Av(V):
    I = np.array([1,1,1])
    i = 1
    j = 1
    while j+3< V.shape[0]:
        while i+3< V.shape[1]:
                R = V[i-1:i+2, j-1:j+2]
                IR = np.matmul(np.matmul(R,I),I)
                R = np.ones_like(R)*IR
        i+=3
    j+=3
    return V

T = np.ones((9,9))
print(Av(T))










def  add(x,y):
     vadd = np.vectorize(add)
     return x+y
R = np.array([1,2,3])
print(vadd(5, R))

#lol, bara använd numpy additio...
u = R +5
print(u)






Q = np.array([[0,-3,5,-7],
              [3,0,13,-17]])
Q[Q > 0] = 0
print(Q)

E = np.array(Q < 0) #E är en boolean array
print(E)
r= Q[E]
print(r)
#or you could use the more elegant syntax:

B = np.array([[True, True,  True,  False],
             [False, False, True, False]])
A = np.zeros((2,4))
AA = []

for i in range(2):
    for j in range(4):
        if Q[i][j] < 0:
            A[i][j] = Q[i][j]
        else:
            A[i][j] = False
print(A)



print(Q.shape)
print(list(Q.shape))









a = np.arange(10)
b = a[[3,7]]
print(b)
c = a[1:3]
print(c)











































































def E(G):
    e = np.zeros((4,4,4,4))
    e[0][1][2][3] = 1
    e[3][0][1][2] = 1
    e[2][3][0][1] = 1
    e[1][2][3][0] = 1
    e[3][2][1][0] = -1
    e[0][3][2][1] = -1
    e[1][0][3][2] = -1
    e[2][1][0][3] = -1
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    if i==j or i==k or  i==l or j==k or j==l or k==l:
                        e[i][j][k][l] = 0
    return np.sqrt(G)*e

print(E(1))
G = np.zeros((4,4))

F = np.array([[0,3,5,7],
              [3,0,13,17],
              [5,-13,0,27],
              [-11,17,27,0]])
A =[]
for k in range(4):
    for l in range(4):
       G[i,j] = np.array(A.append(E(1)[i,j,k,l]*F[k,l]))
print(G)








F = np.array([[0,3,5,7],
              [3,0,13,17],
              [5,-13,0,27],
              [-11,17,27,0]])

































































































