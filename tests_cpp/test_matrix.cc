#include <matrix.h>

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

  //
  //  // test Rows
  // std::cout << "x(1:3,:) = " << x.Rows(4, 7) << std::endl;
  //
  //  //// test Col
  //  std::cout << "x(:,1) = " << x.Col(1) << std::endl;
  //  //
  //  //// test Cols
  //  std::cout << "x(:,1:3) = " << x.Cols(3, 6) << std::endl;
  //
  //  // test Transpose
  //  std::cout << "x^T = " << bla::Transpose(x) << std::endl;

  // test Diag
  std::cout << "x.Diag() = " << x.Diag() << std::endl;

  // set the diagonal to -1
  x.Diag() = -1;
  std::cout << "x = " << x << std::endl;
}