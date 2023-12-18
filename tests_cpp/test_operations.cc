
/*
// #include <lapack_interface.h>
// #include <matrix.h>
#include <cassert>
#include <complex>
#include <iostream>

#include "matrix.h"
#include "vector.h"

namespace bla = Tombino_bla;
// I want to use bla::Vector and bla::Matrix instead of Tombino_bla::Vector and
// Tombino_bla::Matrix
using bla::Matrix;
using bla::ORDERING;
using bla::Vector;

typedef std::complex<double> dcomplex;

template <typename TA, typename TB>
void TestVecVecSum(const Vector<TA>& x, const Vector<TB>& y, bool print = false)
{ typedef decltype(x(0) + y(0)) TRES; Vector<TRES> z(x.Size()); z = x + y; if
(print == true) { std::cout << "Vector Vector Sum " << std::endl; std::cout <<
"x = " << std::endl; std::cout << x << std::endl; std::cout << "y = " <<
std::endl; std::cout << y << std::endl; std::cout << "z = " << std::endl;
    std::cout << z << std::endl;
  }

  for (size_t i = 0; i < x.Size(); i++) {
    TA a = x(i);
    TB b = y(i);
    TRES c = z(i);
    assert(std::abs(a + b - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB>
void TestVecVecSub(const Vector<TA>& x, const Vector<TB>& y, bool print = false)
{ typedef decltype(x(0) - y(0)) TRES; Vector<TRES> z(x.Size()); z = x - y; if
(print == true) { std::cout << "Vector Vector Sub " << std::endl; std::cout <<
"x = " << std::endl; std::cout << x << std::endl; std::cout << "y = " <<
std::endl; std::cout << y << std::endl; std::cout << "z = " << std::endl;
    std::cout << z << std::endl;
  }

  for (size_t i = 0; i < x.Size(); i++) {
    TA a = x(i);
    TB b = y(i);
    TRES c = z(i);
    assert(std::abs(a - b - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecProd(const SCLR& alpha, const Vector<T>& x, bool print = false)
{ typedef decltype(alpha * x(0)) TRES; Vector<TRES> z(x.Size()); z = alpha * x;
  if (print == true) {
    std::cout << "Scalar Vector Mult " << std::endl;
    std::cout << "alpha = " << std::endl;
    std::cout << alpha << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "z = " << std::endl;
    std::cout << z << std::endl;
  }

  for (size_t i = 0; i < x.Size(); i++) {
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha * a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecSum(const SCLR& alpha, const Vector<T>& x, bool print = false) {
  // check that if alpha + x = z using a for loop. If the assert fails, print
the
  // values of x,y,z
  typedef decltype(alpha + x(0)) TRES;
  Vector<TRES> z(x.Size());
  z = alpha + x;
  if (print == true) {
    std::cout << "Scalar Vector Sum " << std::endl;
    std::cout << "alpha = " << std::endl;
    std::cout << alpha << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "z = " << std::endl;
    std::cout << z << std::endl;
  }
  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha + a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecSub(const SCLR& alpha, const Vector<T>& x, bool print = false) {
  // check that if alpha + x = z using a for loop. If the assert fails, print
the
  // values of x,y,z
  typedef decltype(alpha - x(0)) TRES;
  Vector<TRES> z(x.Size());
  z = alpha - x;
  if (print == true) {
    std::cout << "Scalar Vector Sub " << std::endl;
    std::cout << "alpha = " << std::endl;
    std::cout << alpha << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "z = " << std::endl;
    std::cout << z << std::endl;
  }
  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha - a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB>
void TestVecVecOuterProd(const Vector<TA>& x, const Vector<TB>& y, bool print =
false) { typedef decltype(x(0) * y(0)) TRES;
  {
    Matrix<TRES, ORDERING::RowMajor> z(x.Size(), y.Size());

    z = x ^ y;
    if (print == true) {
      std::cout << "Vector Vector Outer Prod: RowMajor " << std::endl;
      std::cout << "x = " << std::endl;
      std::cout << x << std::endl;
      std::cout << "y = " << std::endl;
      std::cout << y << std::endl;
      std::cout << "z = " << std::endl;
      std::cout << z << std::endl;
    }
    for (size_t i = 0; i < x.Size(); i++) {
      for (size_t j = 0; j < y.Size(); j++) {
        // assert with the tolerance of 1e-10
        TA a = x(i);
        TB b = y(j);
        TRES c;
        c = z(i, j);
        assert(std::abs(a * b - c) < 1e-10);
      }
    }
  }

  {
    Matrix<TRES, ORDERING::RowMajor> z(x.Size(), y.Size());
    z = x ^ y;
    if (print == true) {
      std::cout << "Vector Vector Outer Prod: RowMajor " << std::endl;
      std::cout << "x = " << std::endl;
      std::cout << x << std::endl;
      std::cout << "y = " << std::endl;
      std::cout << y << std::endl;
      std::cout << "z = " << std::endl;
      std::cout << z << std::endl;
    }
    for (size_t i = 0; i < x.Size(); i++) {
      for (size_t j = 0; j < y.Size(); j++) {
        // assert with the tolerance of 1e-10
        TA a = x(i);
        TB b = y(j);
        TRES c = z(i, j);
        assert(std::abs(a * b - c) < 1e-10);
      }
    }
  }
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB>
void TestVecVecInnerProd(const Vector<TA>& x, const Vector<TB>& y, bool print =
false) {
  // check that if alpha + x = z using a for loop. If the assert fails, print
the
  // values of x,y,z
  typedef decltype(x(0) * y(0)) TRES;
  TRES c = x * y;
  TRES sum = x(0) * y(0);
  for (size_t i = 0; i < x.Size(); i++) {
    sum = sum + x(i) * y(i);
  }
  if (print == true) {
    std::cout << "Vector Vector Inner Prod " << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "y = " << std::endl;
    std::cout << y << std::endl;
    std::cout << "c = " << std::endl;
    std::cout << c << std::endl;
    std::cout << "sum = " << std::endl;
    std::cout << sum << std::endl;
  }
  assert(std::abs(sum - c) < 1e-10);
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
void TestMatMatSum(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y, bool
print = false) { typedef decltype(X(0, 0) + Y(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), X.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = X + Y;
  Z_col = X + Y;
  if (print == false) {
    std::cout << "Matrix Matrix Sum " << std::endl;
    std::cout << "X = " << std::endl;
    std::cout << X << std::endl;
    std::cout << "Y = " << std::endl;
    std::cout << Y << std::endl;
    std::cout << "Z_row = " << std::endl;
    std::cout << Z_row << std::endl;
    std::cout << "Z_col = " << std::endl;
    std::cout << Z_col << std::endl;
  }

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TB b = Y(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(a + b - c_row) < 1e-10);
      assert(std::abs(a + b - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 9. X - Y
template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
void TestMatMatSub(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y, bool
print = false) { typedef decltype(X(0, 0) - Y(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), X.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = X - Y;
  Z_col = X - Y;
  if (print == true) {
    std::cout << "Matrix Matrix Sub " << std::endl;
    std::cout << "X = " << std::endl;
    std::cout << X << std::endl;
    std::cout << "Y = " << std::endl;
    std::cout << Y << std::endl;
    std::cout << "Z_row = " << std::endl;
    std::cout << Z_row << std::endl;
    std::cout << "Z_col = " << std::endl;
    std::cout << Z_col << std::endl;
  }

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TB b = Y(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(a - b - c_row) < 1e-10);
      assert(std::abs(a - b - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 10. alpha * X
template <typename SCLR, typename TA, ORDERING ORDA>
void TestScalMatProd(const SCLR& alpha, const Matrix<TA, ORDA>& X, bool print =
false) { typedef decltype(alpha * X(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), X.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = alpha * X;
  Z_col = alpha * X;

  if (print == true) {
    std::cout << "Scalar Matrix Mult " << std::endl;
    std::cout << "X = " << std::endl;
    std::cout << X << std::endl;
    std::cout << "Z_row = " << std::endl;
    std::cout << Z_row << std::endl;
    std::cout << "Z_col = " << std::endl;
    std::cout << Z_col << std::endl;
  }

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(alpha * a - c_row) < 1e-10);
      assert(std::abs(alpha * a - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 11. alpha + X
template <typename SCLR, typename TA, ORDERING ORDA>
void TestScalMatSum(const SCLR& alpha, const Matrix<TA, ORDA>& X, bool print =
false) { typedef decltype(alpha + X(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), X.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = alpha + X;
  Z_col = alpha + X;

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(alpha + a - c_row) < 1e-10);
      assert(std::abs(alpha + a - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 12. alpha - X
template <typename SCLR, typename TA, ORDERING ORDA>
void TestScalMatSub(const SCLR& alpha, const Matrix<TA, ORDA>& X, bool print =
false) { typedef decltype(alpha - X(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), X.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = alpha - X;
  Z_col = alpha - X;

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      // print alpha, a, c_row, c_col
      std::cout << "alpha = " << alpha << std::endl;
      std::cout << "a = " << a << std::endl;
      std::cout << "c_row = " << c_row << std::endl;
      std::cout << "c_col = " << c_col << std::endl;
      assert(std::abs(alpha - a - c_row) < 1e-10);
      assert(std::abs(alpha - a - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 13. X * Y
template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
void TestMatMatProd(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y, bool
print = false) { typedef decltype(X(0, 0) * Y(0, 0)) TRES; Matrix<TRES,
ORDERING::RowMajor> Z_row(X.SizeRows(), Y.SizeCols()); Matrix<TRES,
ORDERING::RowMajor> Z_col(X.SizeRows(), Y.SizeCols());

  Z_row = X * Y;
  Z_col = X * Y;

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < Y.SizeCols(); j++) {
      TRES sum = 0;
      for (size_t k = 0; k < X.SizeCols(); k++) {
        sum = sum + X(i, k) * Y(k, j);
      }
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(sum - c_row) < 1e-10);
      assert(std::abs(sum - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 14. X * x
template <typename TA, typename TV, ORDERING ORDA>
void TestMatVecProd(const Matrix<TA, ORDA>& X, const Vector<TV>& x, bool print =
false) { typedef decltype(X(0, 0) * x(0)) TRES; Vector<TRES> z(X.SizeRows()); z
= X * x;

  for (size_t i = 0; i < x.Size(); i++) {
    TRES sum_row_i;
    sum_row_i = X(i, 0) * x(0);
    for (size_t j = 1; j < X.SizeCols(); j++) {
      sum_row_i = sum_row_i + X(i, j) * x(j);
    }
    TRES c_i = z(i);
    assert(std::abs(sum_row_i - c_i) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

// 15. x.Slice()
template <typename T>
void TestVecSlice(const Vector<T>& x, const size_t& first, const size_t& slice,
bool print = false) { Vector<T> y(x.Slice(first, slice)); if (print == true) {
    std::cout << "Vector Slice " << std::endl;
    std::cout << "first = " << std::endl;
    std::cout << first << std::endl;
    std::cout << "slice = " << std::endl;
    std::cout << slice << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "y = " << std::endl;
    std::cout << y << std::endl;
    for (size_t i = 0; i < (x.Size() - first) / slice; i++) {
      // pinrt
      std::cout << "i = " << i << std::endl;
      std::cout << "i * slice + first = " << i * slice + first << std::endl;
      std::cout << "x(i * slice + first) = " << x(i * slice + first) <<
std::endl; std::cout << "y(i) = " << y(i) << std::endl;

      assert(x(i * slice + first) == y(i));
    }
    std::cout << "passed" << std::endl;
  }
}

// 16. x.Range()
template <typename T>
void TestVecRange(const Vector<T>& x, const size_t& first, const size_t& next,
bool print = false) { Vector<T> y(x.Range(first, next)); if (print == true) {
    std::cout << "Vector Range " << std::endl;
    std::cout << "first = " << std::endl;
    std::cout << first << std::endl;
    std::cout << "next = " << std::endl;
    std::cout << next << std::endl;
    std::cout << "x = " << std::endl;
    std::cout << x << std::endl;
    std::cout << "y = " << std::endl;
    std::cout << y << std::endl;
  }
  for (size_t i = 0; i < next - first; i++) {
    assert(x(i + first) == y(i));
  }
  std::cout << "passed" << std::endl;
}

// 17. X.T
template <typename T, ORDERING ORD>
void TestMatTran(const Matrix<T, ORD> X, bool print = false) {
  Matrix<T, ORD> Y_r(Transpose(X));
  Matrix<T, ORD> Y_c(Transpose(X));
  if (print == true) {
    std::cout << "Matrix Transpose " << std::endl;
    std::cout << "X = " << std::endl;
    std::cout << X << std::endl;
    std::cout << "Y_r = " << std::endl;
    std::cout << Y_r << std::endl;
    std::cout << "Y_c = " << std::endl;
  }

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      assert(abs(X(i, j) - Y_r(j, i)) < 1e-10);
      assert(abs(X(i, j) - Y_c(j, i)) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 18. X.Row()
template <typename T, ORDERING ORD>
void TestMatRow(const Matrix<T, ORD> X, const size_t& first, bool print = false)
{ Vector<T> y(X.Row(first)); for (size_t i = 0; i < X.SizeCols(); i++) {
    assert(X(first, i) == y(i));
  }
  std::cout << "passed" << std::endl;
}

// 19. X.Col()
template <typename T, ORDERING ORD>
void TestMatCol(const Matrix<T, ORD> X, const size_t& first, bool print = false)
{ Vector<T> y(X.Col(first)); for (size_t i = 0; i < X.SizeRows(); i++) {
    assert(X(i, first) == y(i));
  }
  std::cout << "passed" << std::endl;
}

// 20. X.Rows()
template <typename T, ORDERING ORD>
void TestMatRows(const Matrix<T, ORD> X, const size_t& first, const size_t&
next, bool print = false) { std::cout << "\nX" << std::endl; std::cout << X <<
std::endl;

  std::cout << "\nX.Rows(first, next)" << std::endl;
  std::cout << X.Rows(first, next) << std::endl;

  Matrix<T, ORDERING::RowMajor> Y(X.Rows(first, next));
  std::cout << "\nY" << std::endl;
  std::cout << Y << std::endl;

  for (size_t i = 0; i < next - first; i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      assert(X(i + first, j) == Y(i, j));
    }
  }
  std::cout << "passed" << std::endl;
}

// 21. X.Cols()
template <typename T, ORDERING ORD>
void TestMatCols(const Matrix<T, ORD> X, const size_t& first, const size_t&
next, bool print = false) { Matrix<T, ORD> Y(X.Cols(first, next)); for (size_t i
= 0; i < X.SizeRows(); i++) { for (size_t j = 0; j < next - first; j++) {
      assert(X(i, j + first) == Y(i, j));
    }
  }
  std::cout << "passed" << std::endl;
}

// 22. Inverse
void TestMatInv(size_t n, bool print = false) {
  // create the matrix with -2 on the diagonal and 1 on the off-diagonal
  Matrix<dcomplex, ORDERING::ColMajor> X(n, n);
  Matrix<dcomplex, ORDERING::ColMajor> I(n, n);
  X = 0;
  I = 0;
  for (size_t i = 0; i < n; i++) {
    X(i, i) = -2;
    I(i, i) = 1;
    if (i < n - 1) {
      X(i, i + 1) = dcomplex(0, 1);
      X(i + 1, i) = dcomplex(0, 1);
    }
  }
  std::cout << "\nX = " << std::endl;
  std::cout << X << std::endl;

  std::cout << "\nI = " << std::endl;
  std::cout << I << std::endl;

  // create the identity matrix
  Matrix<dcomplex, ORDERING::ColMajor> Xinv(n, n);
  Xinv = Inverse(X);

  std::cout << "\nXinv = " << std::endl;
  std::cout << Xinv << std::endl;

  Matrix<dcomplex, ORDERING::ColMajor> Icheck(n, n);
  Icheck = X * Xinv;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      std::cout << (std::abs(I(i, j) - Icheck(i, j))) << std::endl;
      assert(std::abs(I(i, j) - Icheck(i, j)) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 23. Operator()
void TestMatOp(size_t m, size_t n, bool print = false) {
  // here I need to create a matrix of which I know the entries and I need to
check
  // if the entries change in the way that I desired to build the previous tests
  Matrix<double, ORDERING::ColMajor> X(m, n);
  // I fill up the matrix with the following rule:
  // X(i,j) = i*n+j
  // I check that the entries are correct
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) X(i, j) = i * n + j;
  }
  // I check that the entries are correct to do so I print the matrix
  std::cout << "\nX = " << std::endl;
  std::cout << X << std::endl;

  // Now I transpose the matrix and I check that the entries are correct
  std::cout << "\nX.T = " << std::endl;
  std::cout << Transpose(X) << std::endl;

  // Now I change the entries of the matrix and I check that the entries are
correct
  // the last row -1 of the matrix is filled with -10
  for (size_t j = 0; j < n; j++) X(n - 1, j) = -10;
  std::cout << "\nX(n-1,  :) = -10 " << std::endl;
  std::cout << X << std::endl;

  // To check that the above is the same of Rows(n-1,n) I change the entries to
-20 X.Rows(n - 3, n - 1) = -20; std::cout << "\nX.Rows(n-1,n) = -20 " <<
std::endl; std::cout << X << std::endl;

  // now do the same for the columns
  // the last column -1 of the matrix is filled with -30
  std::cout << "\n  m = " << m << std::endl;
  for (size_t i = 0; i < m; i++) X(i, m - 2) = -30;
  std::cout << "\nX(:, n-1) = -30 " << std::endl;
  std::cout << X << std::endl;
  // I check that the entries are correct

  //// test RowSlice()
  // std::cout << "\nX.RowSlice(1, 0 , 3) " << std::endl;
  // std::cout << X.RowSlice(1, 0, 2) << std::endl;
  //
  //// test ColSlice()
  // std::cout << "\nX.ColSlice(2, 1 , 2) " << std::endl;
  // std::cout << X.ColSlice(2, 1, 2) << std::endl;

  // test RowsSlice()
  std::cout << "\nX.RowsSlice(1, 2) " << std::endl;
  std::cout << X.RSlice(1, 2) << std::endl;

  // test ColsSlice()
  std::cout << "\nX.ColsSlice(2, 3) " << std::endl;
  std::cout << X.CSlice(2, 3) << std::endl;

  std::cout << "\nX.RowsSlice(1, 2).CSlice(2, 3) " << std::endl;
  std::cout << X.RSlice(1, 2).CSlice(2, 3) << std::endl;

  // set the matrix to zero

  // std::cout << "\nX.RowsSlice(1,2).CSlice(2,3).Diag() " << std::endl;
  // std::cout << X.RSlice(1, 2).CSlice(2, 3).Diag() << std::endl;
}

void TestDiag(int m, int n)
{
  // create a matrix m x n and initialize with the usual indexing
  Matrix<double, ORDERING::ColMajor> X(m, n);
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++) X(i, j) = i * n + j;
  }
  // print the matrix and diag
  std::cout << "\nX = " << std::endl;
  std::cout << X << std::endl;
  std::cout << "\nX.Diag() = " << std::endl;
  std::cout << X.Diag() << std::endl;

  // extract the diagonal and change it with -2
  X = 0;

  X.Diag() = -2;
  std::cout << "\nX.Diag() = -2 " << std::endl;
  std::cout << X << std::endl;

  X.Diag().Range(1, 3) = -3;
  std::cout << "\nX.Diag().Range(1, 3) = -3 " << std::endl;
  std::cout << X << std::endl;
}

int main() {
  size_t m = 8;
  size_t n = 8;

  size_t l = 8;
  size_t first = 1;  // first <=   n
  size_t slice = 3;  // first <=  slice <= n
  size_t next = 4;   // first <=  next <= n
  bool print = true;

  double alphaR = 2;
  dcomplex alphaC(2, 3);

  dcomplex I(0, 1);

  Vector<double> xR(l), yR(l), zR(l);
  Vector<dcomplex> xC(l), yC(l), zC(l);

  Matrix<double, ORDERING::RowMajor> XR_row(m, n), YR_row(n, l), ZR_row(m, l);
  Matrix<double, ORDERING::RowMajor> XR_col(m, n), YR_col(n, l), ZR_col(m, l);
  Matrix<dcomplex, ORDERING::RowMajor> XC_row(m, n), YC_row(n, l), ZC_row(m, l);
  Matrix<dcomplex, ORDERING::RowMajor> XC_col(m, n), YC_col(n, l), ZC_col(m, l);

  for (size_t i = 0; i < l; i++) {
    xR(i) = i;
    yR(i) = 1;
    xC(i) = dcomplex(i, 1.0);
    yC(i) = dcomplex(1.0, i);
  }

  // fill up the matrices with 1 or 1+I uin the lower triangular part
  for (size_t i = 0; i < m; i++) {
    XR_row(i, i) = 1.0;
    XR_col(i, i) = 1.0;
    XC_row(i, i) = 1.0;
    XC_col(i, i) = 1.0;
    for (size_t j = 0; j < i; j++) {
      XR_row(i, j) = 1.0;
      XR_col(j, i) = 1.0;
      XC_row(i, j) = (1.0 + I);
      XC_col(j, i) = (1.0 + I);
    }
  }
  print = false;
  // 1. x+y
  std::cout << "Test 1: x+y" << std::endl;
  TestVecVecSum(xR, yR, print);
  TestVecVecSum(xC, yC, print);
  TestVecVecSum(xR, yC, print);
  TestVecVecSum(xC, yR, print);

  // 2. x-y
  std::cout << "Test 2: x-y" << std::endl;
  TestVecVecSub(xR, yR, print);
  TestVecVecSub(xC, yC, print);
  TestVecVecSub(xR, yC, print);
  TestVecVecSub(xC, yR, print);

  // 3. alpha*x
  std::cout << "Test 3: alpha*x" << std::endl;
  TestScalVecProd(alphaR, xR, print);
  TestScalVecProd(alphaC, xC, print);
  TestScalVecProd(alphaR, xC, print);
  TestScalVecProd(alphaC, xR, print);

  // 4. alpha + x
  std::cout << "Test 4: alpha + x" << std::endl;
  TestScalVecSum(alphaR, xR, print);
  TestScalVecSum(alphaC, xC, print);
  TestScalVecSum(alphaR, xC, print);
  TestScalVecSum(alphaC, xR, print);
  // 5. alpha - x
  std::cout << "Test 5: alpha - x" << std::endl;
  TestScalVecSub(alphaR, xR, print);
  TestScalVecSub(alphaC, xC, print);
  TestScalVecSub(alphaR, xC, print);
  TestScalVecSub(alphaC, xR, print);

  // 6. x ^ y
  std::cout << "Test 6: x ^ y" << std::endl;
  TestVecVecOuterProd(xR, yR, print);
  TestVecVecOuterProd(xC, yC, print);
  TestVecVecOuterProd(xR, yC, print);
  TestVecVecOuterProd(xC, yR, print);

  // 7. x * y
  std::cout << "Test 7: x * y" << std::endl;
  TestVecVecInnerProd(xR, yR, print);
  TestVecVecInnerProd(xC, yC, print);
  TestVecVecInnerProd(xR, yC, print);
  TestVecVecInnerProd(xC, yR, print);

  // 8. X + Y
  std::cout << "Test 8: X + Y" << std::endl;
  TestMatMatSum(XR_row, YR_row, print);
  TestMatMatSum(XR_col, YR_col, print);
  TestMatMatSum(XC_row, YC_row, print);
  TestMatMatSum(XC_col, YC_col, print);

  // 9. X - Y
  std::cout << "Test 9: X - Y" << std::endl;
  TestMatMatSub(XR_row, YR_row, print);
  TestMatMatSub(XR_col, YR_col, print);
  TestMatMatSub(XC_row, YC_row, print);
  TestMatMatSub(XC_col, YC_col, print);

  // 10. alpha * X
  std::cout << "Test 10: alpha * X" << std::endl;
  TestScalMatProd(alphaR, XR_row, print);
  TestScalMatProd(alphaR, XR_col, print);
  TestScalMatProd(alphaC, XC_row, print);
  TestScalMatProd(alphaC, XC_col, print);

  // 11. alpha + X
  std::cout << "Test 11: alpha + X" << std::endl;
  TestScalMatSum(alphaR, XR_row, print);
  TestScalMatSum(alphaR, XR_col, print);
  TestScalMatSum(alphaC, XC_row, print);
  TestScalMatSum(alphaC, XC_col, print);

  // 12. alpha - X
  std::cout << "Test 12: alpha - X" << std::endl;
  TestScalMatSub(alphaR, XR_row, print);
  TestScalMatSub(alphaR, XR_col, print);
  TestScalMatSub(alphaC, XC_row, print);
  TestScalMatSub(alphaC, XC_col, print);

  // 13. X * Y
  std::cout << "Test 13: X * Y" << std::endl;
  TestMatMatProd(XR_row, YR_row, print);
  TestMatMatProd(XR_col, YR_col, print);
  TestMatMatProd(XC_row, YC_row, print);
  TestMatMatProd(XC_col, YC_col, print);

  // 14. X * x
  std::cout << "Test 14: X * x" << std::endl;
  TestMatVecProd(XR_row, xR, print);
  TestMatVecProd(XR_col, xR, print);
  TestMatVecProd(XC_row, xC, print);
  TestMatVecProd(XC_col, xC, print);
  TestMatVecProd(XR_row, xC, print);
  TestMatVecProd(XR_col, xC, print);
  TestMatVecProd(XC_row, xR, print);
  TestMatVecProd(XC_col, xR, print);

  // test for VectorView
  print = true;
  // 15. X.Slice()
  std::cout << "Test 16: X.Slice()" << std::endl;
  TestVecSlice(xR, first, slice, print);
  TestVecSlice(xC, first, slice, print);

  // 16. X.Range()
  std::cout << "Test 17: X.Range()" << std::endl;
  TestVecRange(xR, first, next, print);
  TestVecRange(xC, first, next, print);

  // Test MatrixView

  // 17. X.Transpose()
  std::cout << "Test 15: X.T" << std::endl;
  TestMatTran(XR_row, print);
  TestMatTran(XR_col, print);
  TestMatTran(XC_row, print);
  TestMatTran(XC_col, print);

  // 18. X.Row()
  std::cout << "Test 18: X.Row()" << std::endl;
  TestMatRow(XR_row, first, print);
  TestMatRow(XR_col, first, print);
  TestMatRow(XC_row, first, print);
  TestMatRow(XC_col, first, print);

  // 19. X.Col()
  std::cout << "Test 19: X.Col()" << std::endl;
  TestMatCol(XR_row, first, print);
  TestMatCol(XR_col, first, print);
  TestMatCol(XC_row, first, print);
  TestMatCol(XC_col, first, print);

  // 20. X.Rows()
  std::cout << "Test 20: X.Rows()" << std::endl;
  TestMatRows(XR_row, first, next, print);
  TestMatRows(XR_col, first, next, print);
  TestMatRows(XC_row, first, next, print);
  TestMatRows(XC_col, first, next, print);

  // 21. X.Cols()
  std::cout << "Test 21: X.Cols()" << std::endl;
  TestMatCols(XR_row, first, next, print);
  TestMatCols(XR_col, first, next, print);
  TestMatCols(XC_row, first, next, print);
  TestMatCols(XC_col, first, next, print);

  // now the tests that require special attention on the entries

  // 22. Inverse
  std::cout << "Test 22: Inverse" << std::endl;
  TestMatInv(2);
  print = true;
  // 23. Operator()
  std::cout << "Test 23: Operator()" << std::endl;
  TestMatOp(10, 9);

  // 24. Diag()

  std::cout << "Test 24: Diag()" << std::endl;
  TestDiag(5, 5);
  TestDiag(5, 10);
  TestDiag(10, 5);
  return 0;
}
*/