from sympy import *
from numpy import *
from sympy import sqrt, exp
from sympy import oo
from sympy.parsing.sympy_parser import parse_expr
x, y, z = symbols('x y z')
init_printing(use_unicode=True)

str = raw_input("Enter integrable function f = ")
f = parse_expr(str)
n = int(raw_input("Enter N = "))
#I
print 'I:'
integral = float(integrate(f, (x, -1, 1)))
Pn2 = 1
Pn1 = x
Pn = 0
for i in xrange(2, n+1):
    Pn = float(2*i-1)/float(i)*x*Pn1-float(i-1)/float(i)*Pn2
    Pn2 = Pn1
    Pn1 = Pn
if n == 0:
    Pn = Pn2
if n == 1:
    Pn = Pn1
Pn = Pn.expand()
X = nroots(Pn, n=15)
print 'X1, ... , Xn:'
print X
A = zeros(n)
g = Pn.diff()
for i in xrange(n):
    A[i] = float(2)/((1-X[i]**2)*(g.subs(x, X[i]))**2)
print 'A1, ... , An:'
print A
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)
#II
print 'II:'
integral = float(integrate(f/sqrt(1-x**2), (x, -1, 1)))
Tn2 = 1
Tn1 = x
Tn = 0
for i in xrange(2, n+1):
    Tn = 2*x*Tn1-Tn2
    Tn2 = Tn1
    Tn1 = Tn
if n == 0:
    Tn = Tn2
if n == 1:
    Tn = Tn1
Tn = Tn.expand()
X = zeros(n)
for i in xrange(n):
    X[i] = cos(((2*i+1)*pi)/(2*n))
print 'X1, ... , Xn:'
print X
A = zeros(n)
for i in xrange(n):
    A[i] = pi/n
print 'A1, ... , An:'
print A
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)
#III
print 'III:'
integral = float(integrate(f*sqrt(1-x**2), (x, -1, 1)))
Tnn = x*2*Tn1-Tn2
Tnn = Tnn.expand()
Tnn = Tnn.diff()
Un = Tnn/(n+1)
X = zeros(n)
for i in xrange(n):
    X[i] = cos((i+1)*pi/(n+1))
print 'X1, ... , Xn:'
print X
A = zeros(n)
for i in xrange(n):
    A[i] = pi/(n+1)*(sin((i+1)*pi/(n+1)))**2
print 'A1, ... , An:'
print A
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)
#IV
print 'IV:'
integral = float(integrate(f*exp(-x**2), (x, -oo, oo)))
u = exp(-x**2)
for i in xrange(n):
    u = u.diff()
Hn = (-1)**n*exp(x**2)*u
Hn = Hn.expand()
X = nroots(Hn, n=15)
print 'X1, ... , Xn:'
print X
g = Hn.diff()
A = zeros(n)
fact = 1
for i in xrange(n):
    fact *= i+1
for i in xrange(n):
    A[i] = 2**(n+1)*fact*sqrt(pi)/(g.subs(x, X[i]))**2
print 'A1, ... , An:'
print A
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)
#V
print 'V:'
a = int(raw_input("Enter alpha, alpha > -1: "))
integral = float(integrate(f*exp(-x)*x**a, (x, 0, oo)))
u = exp(-x)*x**(a+n)
for i in xrange(n):
    u = u.diff()
Ln = (-1)**n*exp(x)*x**(-a)*u
Ln = Ln.expand()
print Ln
X = nroots(Ln, n=15)
print 'X1, ... , Xn:'
print X
A = zeros(n)
G = float(integrate(exp(-x)*x**(a+n), (x, 0, oo)))
g = Ln.diff()
for i in xrange(n):
    A[i] = fact*G/(X[i]*(g.subs(x, X[i]))**2)
print 'A1, ... , An:'
print A
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)
