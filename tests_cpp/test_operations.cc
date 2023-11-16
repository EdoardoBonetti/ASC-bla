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
void TestVecVecSum(const Vector<TA>& x, const Vector<TB>& y) {
  // check that if x+y = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(x(0) + y(0)) TRES;
  Vector<TRES> z(x.Size());
  z = x + y;

  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    TA a = x(i);
    TB b = y(i);
    TRES c = z(i);
    assert(std::abs(a + b - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB>
void TestVecVecSub(const Vector<TA>& x, const Vector<TB>& y) {
  // check that if x-y = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(x(0) - y(0)) TRES;
  Vector<TRES> z(x.Size());
  z = x - y;

  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    TA a = x(i);
    TB b = y(i);
    TRES c = z(i);
    assert(std::abs(a - b - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecProd(const SCLR& alpha, const Vector<T>& x) {
  // check that if alpha*x = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(alpha * x(0)) TRES;
  Vector<TRES> z(x.Size());
  z = alpha * x;

  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha * a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecSum(const SCLR& alpha, const Vector<T>& x) {
  // check that if alpha + x = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(alpha + x(0)) TRES;
  Vector<TRES> z(x.Size());
  z = alpha + x;

  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha + a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename SCLR, typename T>
void TestScalVecSub(const SCLR& alpha, const Vector<T>& x) {
  // check that if alpha + x = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(alpha - x(0)) TRES;
  Vector<TRES> z(x.Size());
  z = alpha - x;

  for (size_t i = 0; i < x.Size(); i++) {
    // assert with the tolerance of 1e-10
    T a = x(i);
    TRES c = z(i);
    assert(std::abs(alpha - a - c) < 1e-10);
  }
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB>
void TestVecVecOuterProd(const Vector<TA>& x, const Vector<TB>& y) {
  typedef decltype(x(0) * y(0)) TRES;
  {
    Matrix<TRES, ORDERING::RowMajor> z(x.Size(), y.Size());

    z = x ^ y;

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
void TestVecVecInnerProd(const Vector<TA>& x, const Vector<TB>& y) {
  // check that if alpha + x = z using a for loop. If the assert fails, print the
  // values of x,y,z
  typedef decltype(x(0) * y(0)) TRES;
  TRES c = x * y;
  TRES sum = x(0) * y(0);
  for (size_t i = 0; i < x.Size(); i++) {
    sum = sum + x(i) * y(i);
  }
  assert(std::abs(sum - c) < 1e-10);
  std::cout << "passed" << std::endl;
}

template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
void TestMatMatSum(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y) {
  typedef decltype(X(0, 0) + Y(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), X.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = X + Y;
  Z_col = X + Y;

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
void TestMatMatSub(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y) {
  typedef decltype(X(0, 0) - Y(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), X.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = X - Y;
  Z_col = X - Y;

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
void TestScalMatProd(const SCLR& alpha, const Matrix<TA, ORDA>& X) {
  typedef decltype(alpha * X(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), X.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = alpha * X;
  Z_col = alpha * X;

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
void TestScalMatSum(const SCLR& alpha, const Matrix<TA, ORDA>& X) {
  typedef decltype(alpha + X(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), X.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), X.SizeCols());

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
void TestScalMatSub(const SCLR& alpha, const Matrix<TA, ORDA>& X) {
  typedef decltype(alpha - X(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), X.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), X.SizeCols());

  Z_row = alpha - X;
  Z_col = alpha - X;

  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      TA a = X(i, j);
      TRES c_row = Z_row(i, j);
      TRES c_col = Z_col(i, j);
      assert(std::abs(alpha - a - c_row) < 1e-10);
      assert(std::abs(alpha - a - c_col) < 1e-10);
    }
  }
  std::cout << "passed" << std::endl;
}

// 13. X * Y
template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
void TestMatMatProd(const Matrix<TA, ORDA>& X, const Matrix<TB, ORDB>& Y) {
  typedef decltype(X(0, 0) * Y(0, 0)) TRES;
  Matrix<TRES, ORDERING::ColMajor> Z_row(X.SizeRows(), Y.SizeCols());
  Matrix<TRES, ORDERING::ColMajor> Z_col(X.SizeRows(), Y.SizeCols());

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
void TestMatVecProd(const Matrix<TA, ORDA>& X, const Vector<TV>& x) {
  typedef decltype(X(0, 0) * x(0)) TRES;
  Vector<TRES> z(X.SizeRows());
  z = X * x;

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

template <typename T, ORDERING ORD>
void TestMatTran(const Matrix<T, ORD> X) {
  Matrix<T, !ORD> Y(X);
  for (size_t i = 0; i < X.SizeRows(); i++) {
    for (size_t j = 0; j < X.SizeCols(); j++) {
      assert(X(i, j) == Y(j, i));
    }
  }
  std::cout << "passed" << std::endl;
}

int main() {
  size_t n = 5;

  double alphaR = 3;
  dcomplex alphaC(3, 4);

  Vector<double> xR(n), yR(n), zR(n);
  Vector<dcomplex> xC(n), yC(n), zC(n);

  Matrix<double, ORDERING::RowMajor> XR_row(n, n), YR_row(n, n), ZR_row(n, n);
  Matrix<double, ORDERING::ColMajor> XR_col(n, n), YR_col(n, n), ZR_col(n, n);
  Matrix<dcomplex, ORDERING::RowMajor> XC_row(n, n), YC_row(n, n), ZC_row(n, n);
  Matrix<dcomplex, ORDERING::ColMajor> XC_col(n, n), YC_col(n, n), ZC_col(n, n);

  for (size_t i = 0; i < n; i++) {
    xR(i) = i;
    yR(i) = 1;
    xC(i) = dcomplex(i, i);
    yC(i) = dcomplex(1, 1);
  }
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      XR_row(i, j) = i + j;
      YR_row(i, j) = 1;
      XR_col(i, j) = i + j;
      YR_col(i, j) = 1;
      XC_row(i, j) = dcomplex(i + j, i + j);
      YC_row(i, j) = dcomplex(1, 1);
      XC_col(i, j) = dcomplex(i + j, i + j);
      YC_col(i, j) = dcomplex(1, 1);
    }
  }
  // 1. x+y
  std::cout << "Test 1: x+y" << std::endl;
  TestVecVecSum(xR, yR);
  TestVecVecSum(xC, yC);
  TestVecVecSum(xR, yC);
  TestVecVecSum(xC, yR);

  // 2. x-y
  std::cout << "Test 2: x-y" << std::endl;
  TestVecVecSub(xR, yR);
  TestVecVecSub(xC, yC);
  TestVecVecSub(xR, yC);
  TestVecVecSub(xC, yR);

  // 3. alpha*x
  std::cout << "Test 3: alpha*x" << std::endl;
  TestScalVecProd(alphaR, xR);
  TestScalVecProd(alphaC, xC);
  TestScalVecProd(alphaR, xC);
  TestScalVecProd(alphaC, xR);

  // 4. alpha + x
  std::cout << "Test 4: alpha + x" << std::endl;
  TestScalVecSum(alphaR, xR);
  TestScalVecSum(alphaC, xC);
  TestScalVecSum(alphaR, xC);
  TestScalVecSum(alphaC, xR);
  // 5. alpha - x
  std::cout << "Test 5: alpha - x" << std::endl;
  TestScalVecSub(alphaR, xR);
  TestScalVecSub(alphaC, xC);
  TestScalVecSub(alphaR, xC);
  TestScalVecSub(alphaC, xR);

  // 6. x ^ y
  std::cout << "Test 6: x ^ y" << std::endl;
  TestVecVecOuterProd(xR, yR);
  TestVecVecOuterProd(xC, yC);
  TestVecVecOuterProd(xR, yC);
  TestVecVecOuterProd(xC, yR);

  // 7. x * y
  std::cout << "Test 7: x * y" << std::endl;
  TestVecVecInnerProd(xR, yR);
  TestVecVecInnerProd(xC, yC);
  TestVecVecInnerProd(xR, yC);
  TestVecVecInnerProd(xC, yR);

  // 8. X + Y
  std::cout << "Test 8: X + Y" << std::endl;
  TestMatMatSum(XR_row, YR_row);
  TestMatMatSum(XR_col, YR_col);
  TestMatMatSum(XC_row, YC_row);
  TestMatMatSum(XC_col, YC_col);

  // 9. X - Y
  std::cout << "Test 9: X - Y" << std::endl;
  TestMatMatSub(XR_row, YR_row);
  TestMatMatSub(XR_col, YR_col);
  TestMatMatSub(XC_row, YC_row);
  TestMatMatSub(XC_col, YC_col);

  // 10. alpha * X
  std::cout << "Test 10: alpha * X" << std::endl;
  TestScalMatProd(alphaR, XR_row);
  TestScalMatProd(alphaR, XR_col);
  TestScalMatProd(alphaC, XC_row);
  TestScalMatProd(alphaC, XC_col);

  // 11. alpha + X
  std::cout << "Test 11: alpha + X" << std::endl;
  TestScalMatSum(alphaR, XR_row);
  TestScalMatSum(alphaR, XR_col);
  TestScalMatSum(alphaC, XC_row);
  TestScalMatSum(alphaC, XC_col);

  // 12. alpha - X
  std::cout << "Test 12: alpha - X" << std::endl;
  TestScalMatSub(alphaR, XR_row);
  TestScalMatSub(alphaR, XR_col);
  TestScalMatSub(alphaC, XC_row);
  TestScalMatSub(alphaC, XC_col);

  // 13. X * Y
  std::cout << "Test 13: X * Y" << std::endl;
  TestMatMatProd(XR_row, YR_row);
  TestMatMatProd(XR_col, YR_col);
  TestMatMatProd(XC_row, YC_row);
  TestMatMatProd(XC_col, YC_col);

  // 14. X * x
  std::cout << "Test 14: X * x" << std::endl;
  TestMatVecProd(XR_row, xR);
  TestMatVecProd(XR_col, xR);
  TestMatVecProd(XC_row, xC);
  TestMatVecProd(XC_col, xC);
  TestMatVecProd(XR_row, xC);
  TestMatVecProd(XR_col, xC);
  TestMatVecProd(XC_row, xR);
  TestMatVecProd(XC_col, xR);

  // test for vector view and matrix view

  // 15. X.Transpose()
  std::cout << "Test 15: X.T" << std::endl;
  TestMatTran(XR_row);
  TestMatTran(XR_col);
  TestMatTran(XC_row);
  TestMatTran(XC_col);
}
