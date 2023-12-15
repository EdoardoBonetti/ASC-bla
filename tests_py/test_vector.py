"""
Test bla module of TomBino: Vector
The test is done with pytest and we compare the results with numpy.
"""

import pickle

import pytest
import numpy as np

from TomBino.bla import Vector


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
    n = 10
    x_tb = Vector(n)
    x_np = np.zeros(n)

    for i in range(n):
        r = np.random.rand()
        x_tb[i] = r
        x_np[i] = r

    # slice the vector

    x_tb_slice = Vector(4)
    x_np_slice = Vector(4)

    x_tb_slice = x_tb[2:6]

    for i in range(4):
        assert x_tb_slice[i] == x_np_slice[i]
