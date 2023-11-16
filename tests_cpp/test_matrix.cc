#include <matrix.h>

#include <iostream>

namespace bla = Tombino_bla;

int main() {
  size_t n = 10;
  bla::Matrix<double, bla::RowMajor> x(n, n), y(n, n), k(n, n);

  for (size_t i = 0; i < x.SizeCols(); i++) {
    for (size_t j = 0; j < x.SizeRows(); j++) {
      x(i, j) = 10 * i + j;
      y(i, i) = -2;
      if ((i != 0) && (i != n - 1)) {
        y(i, i - 1) = 1;
        y(i - 1, i) = 1;
        y(i, 1 + i) = 1;
        y(i + 1, i) = 1;
      }
    }
  }
  // test Row
  std::cout << "x = " << x << std::endl;

  std::cout << "x(1,:) = " << x.Row(1) << std::endl;

  // test Rows
  std::cout << "x(:,1:3) = " << x.Rows(0, 3) << std::endl;

  //// test Col
  std::cout << "x(:,1) = " << x.Col(1) << std::endl;
  //
  //// test Cols
  std::cout << "x(:,1:3) = " << x.Cols(3, 6) << std::endl;

  // test Transpose
  std::cout << "x^T = " << bla::Transpose(x) << std::endl;
}