#include <matrix.h>
#include <vector.h>

#include <iostream>

namespace bla = Tombino_bla;

int main() {
  size_t m = 5;
  size_t n = 10;
  double* data = new double[m * n];
  data[0] = 1;
  bla::MatrixView<double, bla::RowMajor> x(m, n, data);
  x = 2;
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      x(i, j) = i * n + j;
      std::cout << x(i, j) << " ";
    }
    std::cout << std::endl;
  }
  // test Row
  std::cout << "\n\nx = " << x << std::endl;

  //
  std::cout << "x(0,:) = " << x.Rows(0, 2) << std::endl;
  std::cout << "x(1,:) = " << x.Rows(1, 4) << std::endl;
  std::cout << "x(2,:) = " << x.Rows(2, 5) << std::endl;
  std::cout << "x(3,:) = " << x.Rows(3, 7) << std::endl;
  std::cout << "x(4,:) = " << x.Rows(4, 5) << std::endl;
  std::cout << "x(5,:) = " << x.Rows(5, 9) << std::endl;

  // 24. increment

  x.Diag() = -1;

  std::cout << "X = " << std::endl;
  std::cout << x << std::endl;

  // Test increment for vectors
  bla::Vector<double> v(10), w(10);
  v = 1.0;
  w = 2.0;
  v.Range(0, 5) = -10.0;
  std::cout << "v = " << v << std::endl;
  v += w;
  std::cout << "v+= v" << std::endl;
  std::cout << v << std::endl;

  // increment by scalar
  std::cout << "v+= 1.1" << std::endl;
  v *= 1.1;
  std::cout << v << std::endl;

  // Test increment for matrices
  bla::Matrix<double> A(5, 5), B(5, 5);
  A = 1.0;
  B = 2.0;
  A.Rows(2, 4) = -10.0;
  std::cout << "A = " << A << std::endl;
  A += B;
  std::cout << "A+= A" << std::endl;
  std::cout << A << std::endl;

  // increment by scalar
  std::cout << "A+= 1.1" << std::endl;

  A.Cols(1, 4) *= 1.1;
  std::cout << A << std::endl;

  // Test increment for matrices
}