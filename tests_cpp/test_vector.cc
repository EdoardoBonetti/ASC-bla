// #include <lapack_interface.h>
// #include <matrix.h>
#include <matrix.h>
#include <vector.h>

#include <complex>
#include <iostream>

using namespace Tombino_bla;
// using namespace std;

int main() {
  // creates three vectors of size 5
  Vector<double> x(5), y(5), z(5);
  // initialize the vector x with the values 0,1,2,3,4
  for (size_t i = 0; i < 5; i++) {
    x(i) = i;
    y(i) = 1;
  }
  std::cout << "x = " << x << std::endl;
  std::cout << "y = " << y << std::endl;
  std::cout << "z = 3 * x + 3 * y" << std::endl;
  z = 3 * x + 3 * y;
  std::cout << z << std::endl << std::endl;

  typedef std::complex<double> dcomplex;
  dcomplex w;

  dcomplex i1(0, 1);
  w = x * (z + i1 * y);
  std ::cout << w << std::endl;
  std::cout << (x * (z + i1 * y)) << std::endl;
  // initialize the vector x with the values 0,1,2,3,4
}
