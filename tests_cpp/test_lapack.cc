#include <lapack_interface.h>
#include <matrix.h>
#include <vector.h>

#include <iostream>

using std::cout;
using std::endl;
using Tombino_bla::Matrix;
using Tombino_bla::Vector;
using Tombino_bla::ORDERING::ColMajor;

int main() {
  Vector<double> x(5);
  Vector<double> y(5);
  Matrix<double, ColMajor> A(5, 5);
  Matrix<double, ColMajor> B(5, 5);
  Matrix<double, ColMajor> C(5, 5);

  for (size_t i = 0; i < x.Size(); i++) {
    x(i) = i;
    y(i) = 2;
    A(i, i) = i;
    B(i, i) = i;
  }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;

  AddVectorLapack(2, x, y);
  cout << "y+2*x = " << y << endl;

  MultMatMatLapack(A, B, C);
  cout << "A*B = " << C << endl;

  // C = A * B | Lapack;
  //  cout << 2*x+y  Lapack << endl;
}
