import pickle

import pytest
import numpy as np

from TomBino.bla import Matrix, Vector, InnerProduct

def test_matrix_init():
    m, n =  10, 5
    x = Matrix(m, n)
    #assert len(x) == n
    

def test_matrix_set():
    m, n =  10, 5
    x_tb = Matrix(m, n)
    x_np = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            x_tb[i, j] = i + j
            x_np[i, j] = i + j

    for i in range(m):
        for j in range(n):
            assert x_tb[i, j] == x_np[i, j]
            assert x_tb[-i, -j] == x_np[-i, -j]
            assert x_tb[i, -j] == x_np[i, -j]
            assert x_tb[-i, j] == x_np[-i, j]

    
