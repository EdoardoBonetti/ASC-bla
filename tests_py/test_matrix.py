import pickle

import pytest
import numpy as np

from TomBino.bla import Matrix, Vector, InnerProduct


def test_matrix_init():
    m, n = 10, 5
    x = Matrix(m, n)
    assert len(x) == m


def test_matrix_set():
    m, n = 10, 5
    x_tb = Matrix(m, n)
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            x_np[i, j] = i * m + j

    for i in range(m):
        for j in range(n):
            assert x_tb[i, j] == x_np[i, j]
            assert x_tb[-i, -j] == x_np[-i, -j]
            assert x_tb[i, -j] == x_np[i, -j]
            assert x_tb[-i, j] == x_np[-i, j]


def test_matrix_set_slice():
    m, n = 5, 3
    x_tb = Matrix(m, n)
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            x_np[i, j] = i * m + j

    # create s and t slices
    s = slice(0, m, 2)
    t = slice(0, n, 2)

    # i, j must be random integers between -m and m and -n and n
    i = np.random.randint(0, m)
    j = np.random.randint(0, n)
    print(s, t, i, j)
    print("the matrices are:")
    print(x_tb)
    print(x_np)
    print("the column slices are:")
    print(x_tb[:, t])
    print(x_np[:, t])
    # assert np.all(x_tb[s, t] == x_np[s, t])
    # assert np.all(x_tb[i, t] == x_np[i, t])
    # assert np.all(x_tb[s, j] == x_np[s, j])

    print("the row slices are:")
    print(x_tb[s, :])
    print(x_np[s, :])


def test_matrix_add():
    m, n = 10, 5
    x_tb = Matrix(m, n)
    y_tb = Matrix(m, n)

    y_np = np.zeros((m, n))
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            y_tb[i, j] = i * m + j
            x_np[i, j] = i * m + j
            y_np[i, j] = i * m + j

    assert np.all(x_tb + y_tb == x_np + y_np)


def test_matrix_mul():
    m, n = 10, 5
    x_tb = Matrix(m, n)
    y_tb = Matrix(n, m)

    y_np = np.zeros((n, m))
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            y_tb[j, i] = i * m + j
            x_np[i, j] = i * m + j
            y_np[j, i] = i * m + j

    assert np.all(x_tb * y_tb == x_np @ y_np)


def test_matrix_inner_product():
    m, n = 10, 5
    x_tb = Matrix(m, n)
    y_tb = Matrix(m, n)

    y_np = np.zeros((n, m))
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            y_tb[j, i] = i * m + j
            x_np[i, j] = i * m + j
            y_np[j, i] = i * m + j

    print(InnerProduct(x_tb, y_tb))
    sum = 0
    for i in range(m):
        for j in range(n):
            sum += x_tb[i, j] * y_tb[i, j]
    assert InnerProduct(x_tb, y_tb) == sum


def test_matrix_vector_mul():
    m, n = 10, 5
    x_tb = Matrix(m, n)
    y_tb = Vector(n)

    y_np = np.zeros(n)
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * m + j
            x_np[i, j] = i * m + j

    for i in range(n):
        y_tb[i] = i
        y_np[i] = i

    assert np.all(x_tb * y_tb == x_np @ y_np)
    print(x_tb, "\n", y_tb, "\n", x_tb * y_tb)
    print(x_np, "\n", y_np, "\n", x_np @ y_np)


def main():
    test_matrix_inner_product()


if __name__ == "__main__":
    main()
