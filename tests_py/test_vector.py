"""
Test bla module of TomBino: Vector
The test is done with pytest and we compare the results with numpy.
"""

import pickle

import pytest
import numpy as np

from TomBino.bla import Vector


def test_vector_init():
    n = 5
    x = Vector(n)
    assert len(x) == n
    x = 1
    print(x)


def test_vector_set():
    n = 10
    x_tb = Vector(n)
    x_np = np.zeros(n)

    for i in range(n):
        x_tb[i] = i
        x_np[i] = i
    for i in range(n):
        assert x_tb[i] == x_np[i]
        assert x_tb[-i] == x_np[-i]


def test_vector_set_slice():
    n = 10
    x_tb = Vector(n)
    x_np = np.zeros(n)

    for i in range(n):
        x_tb[i] = i
        x_np[i] = i
    for i in range(n):
        assert x_tb[i] == x_np[i]
        assert x_tb[-i] == x_np[-i]

    # create s and t slices
    s = slice(0, n, 2)
    t = slice(0, n, 2)

    # i, j must be random integers between -n and n
    i = np.random.randint(0, n)
    j = np.random.randint(0, n)

    assert np.all(x_tb[s] == x_np[s])
    assert np.all(x_tb[t] == x_np[t])
    assert np.all(x_tb[i] == x_np[i])
    assert np.all(x_tb[j] == x_np[j])


def test_vector_add():
    n = 10
    x_tb = Vector(n)
    y_tb = Vector(n)

    y_np = np.zeros(n)
    x_np = np.zeros(n)

    for i in range(n):
        r = np.random.rand()
        x_tb[i] = r
        x_np[i] = r

        r = np.random.rand()
        y_tb[i] = r
        y_np[i] = r

    z_tb = x_tb + y_tb
    z_np = x_np + y_np

    for i in range(n):
        assert z_tb[i] == z_np[i]

    assert np.all(x_tb + y_tb == x_np + y_np)
    assert np.all(x_tb - y_tb == x_np - y_np)


def main():
    test_vector_init()
    test_vector_set()
    test_vector_add()
    # test_vector_scal_mult()
    # test_vector_slicing()


if __name__ == "__main__":
    main()
