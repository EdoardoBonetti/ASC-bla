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


def main():
    test_matrix_init()
    test_matrix_set()
    test_matrix_set_slice()


if __name__ == "__main__":
    main()
