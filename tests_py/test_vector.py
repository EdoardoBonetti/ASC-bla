"""
Test bla module of TomBino: Vector
The test is done with pytest and we compare the results with numpy.
"""

import pickle

import pytest
import numpy as np

from TomBino.bla import Vector

def test_vector_init():
    n = 10
    x = Vector(n)
    assert len(x) == n

def test_vector_set():
    n = 10
    x_tb = Vector(n)
    print("len(x_tb) = ", len(x_tb))
    print(x_tb)
    x_np = np.zeros(n)
    for i in range(n):
        x_tb[i] = i
        x_np[i] = i
    for i in range(n):
        print(i, x_tb[i], x_np[i])
        assert x_tb[i] == x_np[i]
        print(-i, x_tb[-i], x_np[-i])
        assert x_tb[-i] == x_np[-i]


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


def test_vector_slicing():
    x = Vector(5)
    x[:] = 1
    assert np.array_equal(np.asarray(x), np.ones(5))
    x[1::2] = 2
    assert x[0] == 1
    assert x[3] == 2


def test_vector_add():
    x = Vector(5)
    y = Vector(5)

    for i in range(len(x)):
        x[i] = i
    y[:] = 2

    z = x + y
    assert np.array_equal(np.asarray(z), np.array([2, 3, 4, 5, 6]))


def test_vector_scal_mult():
    x = Vector(5)

    for i in range(len(x)):
        x[i] = i
    z = 2 * x
    assert np.array_equal(np.asarray(z), np.array([0, 2, 4, 6, 8]))