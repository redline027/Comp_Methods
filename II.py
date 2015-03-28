import numpy as np
import sympy
#from sympy import


def Gauss(A, b):
    n, m = A.shape
    l = b.shape[1]
    a = np.zeros((n, n+l))
    a[:, :n] = A
    a[:, n:n+l] = b
    for i in xrange(n):
        b[i, :] /= A[i, i]
        A[i, i+1:m] /= A[i, i]
        for j in xrange(i+1, n):
            b[j, :] -= A[j, i]*b[i, :]
            A[j, i+1:m] -= A[j, i]*A[i, i+1:m]

    c = np.zeros((n, n+l))
    c[:, :n] = A
    c[:, n:n+l] = b

    for i in xrange(n-1, -1, -1):
       for j in xrange(0, i):
           b[j, :] -= A[j, i]*b[i, :]

    return b, a, c


def compGauss(matA, matb):
    n, m = matA.shape
    a = np.zeros((n+1, n+3))
    b = np.zeros((n+1, n+1))
    c = np.zeros((n+1, n+3))
    a[1:n+1, 1:n+1] = matA
    a[1:n+1, n+1] = matb.reshape(3)
    a[1:n+1, n+2] = np.sum(a[1:n+1, 1:n+2], axis=1)

    for k in xrange(1, n+1):
        for s in xrange(n-k+1):
            sum = 0
            for j in xrange(1, k):
                sum += b[k+s, j]*c[j, k]
            b[k+s, k] = a[k+s, k] - sum
        for p in xrange(1, n - k + 3):
            sum = 0
            for j in xrange(1, k, 1):
                   sum += b[k, j]*c[j, k+p]
            c[k, k+p] = (a[k, k+p] - sum)/b[k, k]
        sum = 0
        for j in xrange(k+1, n+2):
            sum += c[k, j]
        c[k, n+2] = 1 + sum

    x = np.zeros(n+1)
    for k in xrange(n):
        sum = 0
        for j in xrange(1, k+1):
            sum += c[n-k, n-j+1]*x[n-j+1]
        x[n-k] = c[n-k, n+1] - sum

    C = np.zeros((n, n+2))
    for i in xrange(n):
        C[i:(n+1), i] = b[(i+1):(n+2), i+1]
        C[i, i+1:n+2] = c[i+1, i+2:n+3]
    A = a[1:n+1, 1:n+3]

    return x[1:4].reshape((3, 1)), A, C


def my_print(x, A, C):
    print "The initial matrix:"
    print A
    print "The final matrix:"
    print C
    print "Solution:"
    print x
    print "Error:"
    n, m = A.shape
    a = A[:, :n]
    b = A[:, n].reshape((n, 1))
    print np.mat(a)*np.mat(x) - np.mat(b)


def det_diag(a):
    n, m = a.shape
    det = 1
    for i in xrange(n):
        det *= a[i, i]
    return det

def norm(A):
    n, m = A.shape
    c = np.zeros(n)
    for i in xrange(n):
        sum = 0
        for k in xrange(n):
            sum += abs(A[i, k])
        c[i] = sum
    return c.max()

def cond_numb(A, A_inv):
    return norm(A) * norm(A_inv)

def II(A, b):
    print "1) Gauss:"
    x, a, c = Gauss(np.copy(A), np.copy(b))
    my_print(x, a, c)
    det = det_diag(c)

    print "2) Comp. Gauss:"
    x, a, c = compGauss(np.copy(A), np.copy(b))
    my_print(x, a, c)

    print "3) Determinant:"
    print det

    print "4) Inverse matrix:"
    A_inv, a, c = Gauss(np.copy(A), np.copy(E))
    print A_inv

    print "5) Solution using the inverse matrix:"
    x = np.mat(A_inv) * np.mat(b)
    print x

    print "6) Condition number:"
    print cond_numb(A, A_inv)

A = np.array([[6.464, 2.180, 3.434],
              [-1.351, 8.474, 5.224],
              [2.489, -0.459, 4.549]])
b = np.array([[54.194],
              [47.793],
              [35.165]])

A1 = np.array([[3.0, 2, -5],
               [2, -1, 3],
               [1, 2, -1]])

A2 = np.array([[1.0, 1, 1],
               [1, 2, 3],
               [1, 4, 5]])

b1 = np.array([[-1.0],
               [13],
               [9]])

E = np.array([[1.0, 0, 0],
              [0, 1.0, 0],
              [0, 0, 1.0]])

II(A, b)




