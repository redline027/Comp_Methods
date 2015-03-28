from sympy import *
from numpy import *
from sympy.parsing.sympy_parser import parse_expr
x, y, z = symbols('x y z')
init_printing(use_unicode=True)

str = raw_input("Enter integrable function f = ")
f = parse_expr(str)

str = raw_input("Enter the weight w = ")
w = parse_expr(str)
a = int(raw_input("Enter the beginning of the interval a = "))
b = int(raw_input("Enter the end of the interval b = "))
n = int(raw_input("Enter N = "))
m = zeros(2*n)
integral = float(integrate(f*w, (x, a, b)))
for i in xrange(2*n):
    m[i] = integrate(w*x**i, (x, a, b))
A = ones((n, n))
for i in xrange(n):
    for j in xrange(n):
        A[i, j] = m[n-1+i-j]
d = ones(n)
for i in xrange(n):
    d[i] = -m[n+i]
c = linalg.solve(A, d)
W = x**n
for i in xrange(n):
    W += c[i]*x**(n-1-i)
kor = zeros((n, 2))
h = float(raw_input("Enter step for bisection, h = "))
cur = a
i = 0
while cur+h <= b:
    if W.subs(x, cur)*W.subs(x, cur+h) <= 0:
        kor[i] = [cur, cur+h]
        i += 1
    cur += h
X = zeros(n)
for i in xrange(n):
    h = kor[i, 1] - kor[i, 0]
    while h > 0.000000002:
        cur = (kor[i, 1] + kor[i, 0])/2
        if W.subs(x, kor[i, 0])*W.subs(x, cur) <= 0:
            kor[i, 1] = cur
        else:
            kor[i, 0] = cur
        h = kor[i, 1] - kor[i, 0]
    X[i] = (kor[i, 1] + kor[i, 0])/2
C = ones((n, n))
for i in xrange(n):
    for j in xrange(n):
        C[i, j] = X[j]**i
d = ones(n)
for i in xrange(n):
    d[i] = m[i]
A = linalg.solve(C, d)
value = 0
for i in xrange(n):
    value += A[i]*f.subs(x, X[i])
if abs(value - m[2*n-1]) < 0.000000001:
    print "Function is", x**(2*n-1)
print "Integral is equal to ", integral
print "The approximate value of the integral is equal to ", value
print "Error is equal to ", abs(integral - value)