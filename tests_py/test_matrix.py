from TomBino.bla import Vector, VectorInt, VectorComplex
from TomBino.bla import Matrix, MatrixComplex

mc = MatrixComplex(3, 3)
mc[0, 0] = 1 + 1j
print(mc)
print(mc[0, 1])

mc2 = MatrixComplex(10, 7)
for I in range(10):
    for J in range(7):
        mc2[I, J] = I + J * 1j

print(mc2)
print(mc2[0, 1])
print(mc2[0, 2])

mc2[1:, 0:4] = 99 + 99j
print(mc2)

quit()


x = Vector(3)
y = Vector(3)

print("x =", x)
print("y =", y)

x[:] = 1
print("x =", x)
n = 10
X = Matrix(n, n)
Y = Matrix(n, n)
for i in range(n):
    for j in range(n):
        X[i, j] = 10 * i + j

print(X)


# print (X)
# print (X.T)
#
# print(X.inv)

X[3, 0:10:4] = 100

print((X))

x = Vector(10)
for i in range(10):
    x[i] = 10
x[0:4] = 3
print(x)

X = -2 * X

# print(X[0:3, 0:3]* X[0, 0:3])

# Measure timings for Matrix-Matrix multiplicaton called from Python
# (width=height=n, with n=10, n=100, n=1000).
# Split times into actual C++ computations, and overhead due to Python wrapping.
# Compare with Numpy timings in milliseconds.
#
import time
import numpy as np

T_cpp = []
T_np = []

for n in [100 * n for n in range(1, 10)]:
    X = Matrix(n, n)
    Y = Matrix(n, n)
    for i in range(n):
        for j in range(n):
            X[i, j] = 10 * i + j
            Y[i, j] = 10 * i + j

    t0 = time.time()
    Z = X * Y
    t1 = time.time()
    print("n = %4d: My C++ time = %6.3f ms" % (n, (t1 - t0)))
    T_cpp.append(t1 - t0)

    Xnp = np.asarray(X)
    Ynp = np.asarray(Y)
    t2 = time.time()
    Znp = np.dot(Xnp, Ynp)
    t3 = time.time()
    print("n = %4d: Numpy time = %6.3f ms" % (n, (t3 - t2)))
    T_np.append(t3 - t2)

# Plot timings log scale

import matplotlib.pyplot as plt

dims = np.log(np.asarray([100 * n for n in range(1, 14)]))
plt.figure()
plt.plot(dims, np.log(T_cpp), label="C++")
plt.plot(dims, np.log(T_np), label="Numpy")
plt.legend()
plt.xlabel("n")
plt.ylabel("time [s]")
plt.title("Matrix-Matrix multiplication")
format = "svg"
plt.show()
# plt.savefig("test_vector.svg",format=format, dpi=1200)
