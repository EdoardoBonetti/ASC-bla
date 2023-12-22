"""
Test bla module of TomBino: Matrix
The test is done with pytest and we compare the results with numpy.
"""

import pickle

# import pytest
import numpy as np

from TomBino.bla import Matrix, Vector, InnerProduct, Inverse, Transpose


def test_matrix_init():
    """Test the initialization of a matrix"""
    m, n = 10, 5
    x = Matrix(m, n)
    assert len(x) == m


def test_matrix_set():
    """Test the setting of a matrix"""
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
    """Test the setting of a matrix slice"""
    m, n = 5, 3
    x_tb = Matrix(m, n)
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i * n + j
            x_np[i, j] = i * n + j

    # create s and t slices
    s = slice(0, m, 2)
    t = slice(0, n, 2)

    # i, j must be random integers between -m and m and -n and n
    i = np.random.randint(0, m)
    j = np.random.randint(0, n)

    assert np.all(np.asarray(x_tb[i, t]) == x_np[i, t])
    assert np.all(np.asarray(x_tb[s, j]) == x_np[s, j])
    # assert np.all(np.asarray(x_tb[s, t]) == x_np[s, t])
    print(x_tb)
    print(x_np)

    print(x_tb[s, t])
    print(x_np[s, t])



def test_matrix_add():
    """Test the addition of two matrices"""
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
    """Test the multiplication of two matrices"""
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
    """Test the inner product of two matrices"""
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

    total = 0
    for i in range(m):
        for j in range(n):
            total += x_tb[i, j] * y_tb[i, j]
    assert InnerProduct(x_tb, y_tb) == total


def test_matrix_vector_mul():
    """Test the multiplication of a matrix and a vector"""
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



# def test_matrix_pickle():
#    """Test the pickle of a matrix"""
#    m, n = 10, 5
#    x_tb = Matrix(m, n)
#    y_tb = Vector(n)
#
#    y_np = np.zeros(n)
#    x_np = np.zeros((m, n))
#
#    for i in range(m):
#        for j in range(n):
#            x_tb[i, j] = i * m + j
#            x_np[i, j] = i * m + j
#
#    for i in range(n):
#        y_tb[i] = i
#        y_np[i] = i
#
#    with open("test_matrix.pickle", "wb") as f:
#        pickle.dump(x_tb, f)
#
#    with open("test_matrix.pickle", "rb") as f:
#        x_tb = pickle.load(f)
#
#    assert np.all(x_tb * y_tb == x_np @ y_np)
#    print(x_tb, "\n", y_tb, "\n", x_tb * y_tb)
#    print(x_np, "\n", y_np, "\n", x_np @ y_np)


def test_matrix_transpose():
    """Test the transpose of a matrix"""
    m, n = 10, 5
    x_tb = Matrix(m, n)
    x_np = np.zeros((m, n))

    for i in range(m):
        for j in range(n):
            r = np.random.rand()
            x_tb[i, j] = r
            x_np[i, j] = r
    assert np.all(Transpose(x_tb) == x_np.T)


def test_matrix_inverse():
    """Test the inverse of a matrix"""
    for n in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
        x_tb = Matrix(n, n)
        x_np = np.zeros((n, n))

        for i in range(n):
            for j in range(n):
                r = np.random.rand()
                x_tb[i, j] = r
                x_np[i, j] = r
        assert np.all(np.linalg.inv(x_np) - np.asarray(Inverse(x_tb)) < 1e-10)


def test_matrix_vector_mul():
    """Test the multiplication of a matrix and a vector"""
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


# def test_matrix_pickle():
#    """Test the pickle of a matrix"""
#    m, n = 10, 5
#    x_tb = Matrix(m, n)
#    y_tb = Vector(n)
#
#    y_np = np.zeros(n)
#    x_np = np.zeros((m, n))
#
#    for i in range(m):
#        for j in range(n):
#            x_tb[i, j] = i * m + j
#            x_np[i, j] = i * m + j
#
#    for i in range(n):
#        y_tb[i] = i
#        y_np[i] = i
#
#    # pickle is made by a tuple
#    with open("test_matrix.pickle", "wb") as f:
#        pickle.dump((x_tb, y_tb), f)
#
#    with open("test_matrix.pickle", "rb") as f:
#        x_tb = pickle.load(f)
#
#    assert np.all(x_tb == x_np)
#    assert np.all(y_tb == y_np)
#

def main():
    """Run all tests"""
    test_matrix_inverse()
    test_matrix_transpose()
    test_matrix_init()
    test_matrix_set()
    test_matrix_set_slice()
    test_matrix_add()
    test_matrix_mul()
    test_matrix_inner_product()
    test_matrix_vector_mul()
    # test_matrix_pickle()
    print("All tests passed!")


if __name__ == "__main__":
    main()
