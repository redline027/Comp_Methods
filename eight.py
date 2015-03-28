from sympy import *
from numpy import *
from sympy import sqrt, exp
from sympy.parsing.sympy_parser import parse_expr
x, y, z = symbols('x y z')
init_printing(use_unicode=True)

print "Equation is y' = y - x**2"
x0 = float(raw_input("Enter x0 = "))
y0 = float(raw_input("Enter y(x0) = "))
h = float(raw_input("Enter h = "))
N = int(raw_input("Enter N = "))

f = y - x**2
f1 = y - 2*x
f2 = y - 2
fn = y

a = float(y0-x0**2-2*x0-2)/exp(x0)
sol = a*exp(x)+x**2+2*x+2

print "Exact solution:"
for i in range(-2, N+1):
    print "k =", i, " x =", x0+h*i, " y =", sol.subs(x, x0+h*i)

#Teaylor
def Taylor(a):
    T = 0 #result
    i = 0 #add ==0
    j = 0 #add !=0
    p = 0 #derivative
    fact = 6
    while True:
        if i == 5 or j == 5:
            break
        if p == 0:
            if y0 != 0:
                T += y0
#               print T
                j += 1
            p += 1
        elif p == 1:
            der = f.subs(x, x0).subs(y, y0)
            add = der*(a - x0)
            if add != 0:
#                print add
                T += add
                j += 1
            p += 1
        elif p == 2:
            der = f1.subs(x, x0).subs(y, der)
            add = der*(a - x0)**2 / 2
            if add != 0:
#                print add
                T += add
                j += 1
            p += 1
        elif p == 3:
            der = f2.subs(y, der)
            add = der*(a-x0)**3 / 6
            if add != 0:
#                print add
                T += add
                j += 1
            p += 1
        else:
            fact *= p
            add = der*(a-x0)**p/fact
            if add != 0:
#                print add
                T += add
                j += 1
            else:
                i += 1
    return T
Y = zeros(5)
print "Taylor:"
for i in range(-2, 3):
    Y[i+2] = Taylor(x0+h*i)
    print "k =", i, "x =", x0+h*i, "y =", Y[i+2], "error =", abs(sol.subs(x, x0+h*i)-Y[i+2])
#print "k =", 1, "x =", x0+h, "y =", Taylor(x0+h), "error =", abs(sol.subs(x, x0+h)-Taylor(x0+h))

#Adams
Ad = zeros([N+3, 7])
for i in range(N+3):
    Ad[i, 0] = x0+h*(i-2)
for i in range(5):
    Ad[i, 1] = Y[i]
#    Ad[i, 1] = sol.subs(x, Ad[i, 0])
    Ad[i, 2] = h*f.subs(x, Ad[i, 0]).subs(y, Ad[i, 1])
for i in range(3, 7):
    for j in range(4, -1+i-2, -1):
        Ad[j, i] = Ad[j, i-1] - Ad[j-1, i-1]
for i in range(5, N+3):
    Ad[i, 1] = Ad[i-1, 1] + Ad[i-1, 2] + Ad[i-1, 3]/2 + Ad[i-1, 4]*5/12 + Ad[i-1, 5]*3/8 + Ad[i-1, 6]*251/720
    Ad[i, 2] = h*f.subs(x, Ad[i, 0]).subs(y, Ad[i, 1])
    for j in range(3, 7):
        Ad[i, j] = Ad[i, j-1] - Ad[i-1, j-1]
print "Adams:"

for i in range(N+3):
    print "k =", i-2, "x =", Ad[i, 0], "y =", Ad[i, 1], "error =", abs(sol.subs(x, Ad[i, 0])-Ad[i, 1])

#Runge
yn = y0
xn = x0
print "Runge-Kutta:"
for i in range(1, N+1):
    k1 = h*f.subs(x, xn).subs(y, yn)
    k2 = h*f.subs(x, xn+h/2).subs(y, yn+k1/2)
    k3 = h*f.subs(x, xn+h/2).subs(y, yn+k2/2)
    k4 = h*f.subs(x, xn+h).subs(y, yn+k3)
    yn += (k1+2*k2+2*k3+k4)*1/6
    xn += h
    print "k =", i, "x =", xn, "y =", yn, "error =", abs(sol.subs(x, xn)-yn)

#Euler
xn = x0
E1 = y0
E2 = y0
E3 = y0
print "Euler:"
for i in range (1, N+1):
    E1 += h*f.subs(x, xn).subs(y, E1)

    Y2 = E2 + f.subs(x, xn).subs(y, E2)*h/2
    E2 += f.subs(x, xn+h/2).subs(y, Y2)*h

    Y3 = E3 + h*f.subs(x, xn).subs(y, E3)
    E3 += (f.subs(x, xn).subs(y, E3)+f.subs(x, xn+h).subs(y, Y3))*h/2
    xn += h
    print "k =", i, "x =", xn
    print "1: y =", E1, "error =", abs(sol.subs(x, xn)-E1)
    print "2: y =", E2, "error =", abs(sol.subs(x, xn)-E2)
    print "3: y =", E3, "error =", abs(sol.subs(x, xn)-E3)