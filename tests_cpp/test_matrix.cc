#include <matrix.h>

#include <iostream>

namespace bla = Tombino_bla;

int main() {
  size_t m = 5;
  size_t n = 10;
  double* data = new double[m * n];
  data[0] = 1;
  bla::MatrixView<double, bla::ColMajor> x(m, n, data);
  x = 2;
  // for (size_t i = 0; i < m; i++) {
  //   for (size_t j = 0; j < n; j++) {
  //     x(i, j) = i * m + j;
  //   }
  // }
  //   // test Row
  std::cout << "\n\nx = " << x << std::endl;
  //
  //  std::cout << "x(1,:) = " << x.Row(1) << std::endl;
  //
  //  // test Rows
  //  std::cout << "x(1:3,:) = " << x.Rows(4, 7) << std::endl;
  //
  //  //// test Col
  //  std::cout << "x(:,1) = " << x.Col(1) << std::endl;
  //  //
  //  //// test Cols
  //  std::cout << "x(:,1:3) = " << x.Cols(3, 6) << std::endl;
  //
  //  // test Transpose
  //  std::cout << "x^T = " << bla::Transpose(x) << std::endl;
}