"""
Test bla module of TomBino: Vector
The test is done with pytest and we compare the results with numpy.
"""

import pickle

import pytest
import numpy as np

from TomBino.bla import Vector, InnerProduct


# define function random_vector
def random_vectors(n: int) -> (Vector, np.ndarray):
    """Return a random vectors of length n"""
    x = Vector(n)
    y = np.zeros(n)
    for i in range(n):
        r = np.random.rand()
        x[i] = r
        y[i] = r
    return x, y


def random_slice(n: int) -> slice:
    """Return a random slice of Slice(a,b,c):
    a < b  and c is random integer between 1 and b-a
    """

    a = 0
    b = 0
    while a == b:
        a = np.random.randint(0, n)
        b = np.random.randint(0, n)
    if a > b:
        a, b = b, a
    if b - a == 1:
        c = 1
    else:
        c = np.random.randint(1, b - a)
    return slice(a, b, c)


class Test_Vector:

    n = np.random.randint(1, 100)
    x_tb, x_np = random_vectors(n)

    def test_init(self):
        assert len(self.x_tb) == self.n

    def test_set_get_component(self):
        """Test the setter and getter of a vector component"""
        for i in range(self.n):
            assert self.x_tb[i] == self.x_np[i]
            assert self.x_tb[-i] == self.x_np[-i]

        for i in range(self.n):
            r = np.random.rand()
            self.x_tb[i] = r
            self.x_np[i] = r

        for i in range(self.n):
            assert self.x_tb[i] == self.x_np[i]
            assert self.x_tb[-i] == self.x_np[-i]

    def test_set_get_slice(self):
        """Test setter and getter of a vector slice"""

        # create s and t slices
        s = random_slice(self.n)
        t = random_slice(self.n)

        # i, j must be random integers between -n and n
        i = np.random.randint(0, self.n)
        j = np.random.randint(0, self.n)

        assert np.all(self.x_tb[i] == self.x_np[i])
        assert np.all(self.x_tb[j] == self.x_np[j])
        assert np.all(np.asarray(self.x_tb[s]) == self.x_np[s])
        assert np.all(self.x_tb[t] == self.x_np[t])

    @pytest.mark.skip(reason="Not implemented yet")
    def test_vector_add():
        n = 20
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

    @pytest.mark.skip(reason="Not implemented yet")
    def test_vector_scal_mult():
        n = 10
        x_tb = Vector(n)
        x_np = np.zeros(n)

        for i in range(n):
            r = np.random.rand()
            x_tb[i] = r
            x_np[i] = r

        r = np.random.rand()
        x_tb *= r
        x_np *= r

        for i in range(n):
            assert x_tb[i] == x_np[i]

    @pytest.mark.skip(reason="Not implemented yet")
    def test_vector_inner_product():
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

        assert x_tb * y_tb == x_np.dot(y_np)
        assert InnerProduct(x_tb, y_tb) == x_np.dot(y_np)


def main():
    Test_Vector.test_init()
    # test_vector_set()
    # test_vector_add()
    # test_vector_scal_mult()
    # test_vector_inner_product()
    print("All tests passed!")

if __name__ == "__main__":
    main()
