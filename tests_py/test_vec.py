"""
Test bla module of TomBino: Vec
The test is done with pytest and we compare the results with numpy.
"""

# import pickle
# import pytest
import numpy as np

from TomBino.bla import Vec3


def test_vec_init():
    """Test the initialization of a vec"""
    x = Vec3()
    assert len(x) == 3
    x = 1
    print(x)


# def test_vec_set():
#    """Test the setting of a vec"""
#
#    n = 10
#    x_tb = Vec3(n)
#    x_np = np.zeros(n)
#
#    for i in range(n):
#        x_tb[i] = i
#        x_np[i] = i
#    for i in range(n):
#        assert x_tb[i] == x_np[i]
#        assert x_tb[-i] == x_np[-i]
#
#
# def test_vec_set_slice():
#    """Test the setting of a vec slice"""
#    n = 10
#    x_tb = Vec3(n)
#    x_np = np.zeros(n)
#
#    for i in range(n):
#        x_tb[i] = i
#        x_np[i] = i
#    for i in range(n):
#        assert x_tb[i] == x_np[i]
#        assert x_tb[-i] == x_np[-i]
#
#    # create s and t slices
#    s = slice(0, n, 2)
#    t = slice(0, n, 2)
#
#    # i, j must be random integers between -n and n
#    i = np.random.randint(0, n)
#    j = np.random.randint(0, n)
#
#    assert np.all(x_tb[s] == x_np[s])
#    assert np.all(x_tb[t] == x_np[t])
#    assert np.all(x_tb[i] == x_np[i])
#    assert np.all(x_tb[j] == x_np[j])
#
#
# def test_vec_add():
#    """Test the addition of two vecs"""
#    n = 10
#    x_tb = Vec3(n)
#    y_tb = Vec3(n)
#
#    y_np = np.zeros(n)
#    x_np = np.zeros(n)
#
#    for i in range(n):
#        r = np.random.rand()
#        x_tb[i] = r
#        x_np[i] = r
#
#        r = np.random.rand()
#        y_tb[i] = r
#        y_np[i] = r
#
#    z_tb = x_tb + y_tb
#    z_np = x_np + y_np
#
#    for i in range(n):
#        assert z_tb[i] == z_np[i]
#
#    assert np.all(x_tb + y_tb == x_np + y_np)
#    assert np.all(x_tb - y_tb == x_np - y_np)
#
#
# def test_vec_scal_mult():
#    """Test the scalar multiplication of a vec"""
#    n = 10
#    x_tb = Vec3(n)
#    x_np = np.zeros(n)
#
#    for i in range(n):
#        r = np.random.rand()
#        x_tb[i] = r
#        x_np[i] = r
#
#    r = np.random.rand()
#    x_tb *= r
#    x_np *= r
#
#    for i in range(n):
#        assert x_tb[i] == x_np[i]
#
#
# def test_vec_inner_product():
#    """Test the inner product of two vecs"""
#    n = 10
#    x_tb = Vec3(n)
#    y_tb = Vec3(n)
#
#    y_np = np.zeros(n)
#    x_np = np.zeros(n)
#
#    for i in range(n):
#        r = np.random.rand()
#        x_tb[i] = r
#        x_np[i] = r
#
#        r = np.random.rand()
#        y_tb[i] = r
#        y_np[i] = r
#
#    assert x_tb * y_tb == x_np.dot(y_np)
#    assert InnerProduct(x_tb, y_tb) == x_np.dot(y_np)

def main():
    """Main function"""
    test_vec_init()


if __name__ == "__main__":
    main()
