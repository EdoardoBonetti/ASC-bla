#include <lapack_interface.h>
#include <matrix.h>
#include <vector.h>

#include <iostream>

using namespace ASC_bla;
using namespace std;

int main() {
  Vector<double> x(5);
  Vector<double> y(5);

  for (int i = 0; i < x.Size(); i++) {
    x(i) = i;
    y(i) = 2;
  }

  cout << "x = " << x << endl;
  cout << "y = " << y << endl;

  AddVectorLapack(2, x, y);
  cout << "y+2*x = " << y << endl;

  namespace bla = ASC_bla;

  Matrix<double, bla::RowMajor> A(5, 5);
  Matrix<double, bla::RowMajor> B(5, 5);
  Matrix<double, bla::ColMajor> C(5, 5);

  for (int i = 0; i < A.SizeCols(); i++) {
    for (int j = 0; j < A.SizeRows(); j++) {
      A(i, j) = i;
      B(i, j) = j;
    }
  }

  cout << "A = " << A << endl;
  cout << "B = " << B << endl;

  MultMatMatLapack(A, B, C);
  cout << "C = A*B = " << C << endl;
}
