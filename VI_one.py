from sympy import *
from numpy import *
from sympy import sqrt, exp, log
from sympy.parsing.sympy_parser import parse_expr
x, y, z = symbols('x y z')
init_printing(use_unicode=True)

def Sol(N, A, B, C, D):
    M = zeros([N+1, N+1])
    for i in range(N):
        M[i, i] = B[i]
        M[i+1, i] = A[i+1]
        M[i, i+1] = C[i]
    M[N, N] = B[N]
    print "Matrix:"
    for i in range(N+1):
        for j in range(N+1):
            print '%f' % M[i, j],
        print D[i]
    m = zeros(N+1)
    k = zeros(N+1)
    m[0] = -C[0]/B[0]
    k[0] = D[0]/B[0]
    for i in range(N):
        m[i+1] = -C[i]/(A[i]*m[i]+B[i])
        k[i+1] = (D[i]-A[i]*k[i])/(A[i]*m[i]+B[i])
    print "Coefficients:"
    print "M:", m
    print "K:", k
    y = zeros(N+1)
    y[N] = (D[N]-k[N]*A[N])/(A[N]*m[N]+B[N])
    for i in range(N-1, -1, -1):
        y[i] = m[i+1]*y[i+1]+k[i+1]
    print "Solution:", y
    print "Error:"
    print M*mat([y]).T-mat(D).T

"""N = 4
A = array([0, -1, 2, 1, 2])
B = array([1, 1, -2, -2, 2])
C = array([1, -1, 1, 1, 0])
D = array([0, -3, -4, 2, 2])
Sol(N, A, B, C, D)"""

a0 = 0.6
a1 = -1
a = 0
A = 0
b0 = 0.4
b1 = 1
b = 1
B = 0
p = 1 + x*0
q = log(x+2)
r = -x
f = x+2

N = int(raw_input("Enter N = "))

h = float(b-a)/N
X = zeros(N+1)
P = zeros(N+1)
Q = zeros(N+1)
R = zeros(N+1)
F = zeros(N+1)
for i in range(N+1):
    X[i] = a+i*h
    P[i] = p.subs(x, X[i])
    Q[i] = q.subs(x, X[i])
    R[i] = r.subs(x, X[i])
    F[i] = f.subs(x, X[i])
#I
MA = zeros(N+1)
MB = zeros(N+1)
MC = zeros(N+1)
MD = zeros(N+1)
MA[0] = 0
MB[0] = h*a0-a1
MC[0] = a1
MD[0] = h*A
MA[N] = -b1
MB[N] = h*b0+b1
MC[N] = 0
MD[N] = h*B
for i in range(1, N):
    MA[i] = P[i]-h*Q[i]/2
    MB[i] = -2*P[i]+R[i]*h*h
    MC[i] = P[i]+h*Q[i]/2
    MD[i] = h*h*F[i]
print "First:"
Sol(N, MA, MB, MC, MD)

#II
MA[0] = 0
MB[0] = 2*h*a0+a1*(MA[1]/MC[1]-3)
MC[0] = a1*(MB[1]/MC[1]+4)
MD[0] = 2*h*A+a1*MD[1]/MC[1]
MA[N] = -b1*(4+MB[N-1]/MA[N-1])
MB[N] = 2*h*b0+b1*(3-MC[N-1]/MA[N-1])
MC[N] = 0
MD[N] = 2*h*B-b1*MD[N-1]/MA[N-1]
print "Second:"
Sol(N, MA, MB, MC, MD)
