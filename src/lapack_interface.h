// Purpose: Interface to LAPACK library.
// Date: 10.11.2023

// To Do:
// - add documentation
// - add tests
// - add complex case

#ifndef FILE_LAPACK_INTERFACE_H
#define FILE_LAPACK_INTERFACE_H

#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "expression.h"
#include "matrix.h"
#include "vector.h"

typedef int integer;
typedef integer logical;
typedef float real;
typedef double doublereal;
typedef std::complex<float> singlecomplex;
typedef std::complex<double> doublecomplex;

// Windows SDK defines VOID in the file WinNT.h
#ifndef VOID
typedef void VOID;
#endif
typedef int ftnlen;
typedef int L_fp;  // ?

extern "C" {
#include <clapack.h>

#include "clapack.h"
}

namespace Tombino_bla {

// add vecs
template <typename SX, typename SY>
void AddVectorLapack(double alpha, VectorView<double, SX> x, VectorView<double, SY> y) {
  integer n = x.Size();
  integer incx = x.Dist();
  integer incy = y.Dist();
  daxpy_(&n, &alpha, &x(0), &incx, &y(0), &incy);
};

// mult vecs
template <ORDERING OA, ORDERING OB>
void MultMatMatLapack(MatrixView<double, OA> a, MatrixView<double, OB> b, MatrixView<double, ColMajor> c) {
  char transa = (OA == ColMajor) ? 'N' : 'T';
  char transb = (OB == ColMajor) ? 'N' : 'T';

  integer n = c.SizeRows();
  integer m = c.SizeCols();
  integer k = a.SizeCols();
  // if (n == 0 || m == 0) return 0;

  double alpha = 1.0;
  double beta = 0;
  integer lda = a.DistCols();
  integer ldb = b.DistCols();
  integer ldc = c.DistCols();

  int errcode = dgemm_(&transa, &transb, &n, &m, &k, &alpha, a.Data(), &lda, b.Data(), &ldb, &beta, c.Data(), &ldc);
  // if (errcode != 0) throw exception("Lapack-dgemm error " + std::to_string(errcode));
}

// mult vecs 2
template <ORDERING OA, ORDERING OB>
void MultMatMatLapack(MatrixView<double, OA> a, MatrixView<double, OB> b, MatrixView<double, RowMajor> c) {
  MultMatMatLapack(Transpose(b), Transpose(a), Transpose(c));
}

class T_Lapack {};
static constexpr T_Lapack Lapack;

// define a class LapackMultExpr containing the factors A and B
template <typename TA, typename TB>
class LapackMultExpr {
  TA A_;
  TB B_;

 public:
  LapackMultExpr(const TA& A, const TB& B) : A_(A), B_(B) {}

  const TA& A() const { return A_; }
  const TB& B() const { return B_; }
};
// operator| acting on :
// - MultExpr<MatrixView,MatrixView>
// - Lapack (T_Lapack)
// The operator| returns an object of type LapackMultExpr which contains the
// factors A and B
//
////template <ORDERING OA, ORDERING OB>
// auto operator|(const ProdMatExpr<MatrixView<double, OA>, MatrixView<double,
// OB>>& expr, T_Lapack) {
//   return LapackMultExpr<MatrixView<double, OA>, MatrixView<double,
//   OB>>(expr.A(), expr.B());
// }
//// define an assignment operator for the class MatrixView taking such an
/// LapackMultExpr type, /     and calling the MultMatMatLapack function
// template <ORDERING OA, ORDERING OB>
// MatrixView<double, ColMajor>& operator=(MatrixView<double, ColMajor>& C,
//                                         const
//                                         LapackMultExpr<MatrixView<double,
//                                         OA>, MatrixView<double, OB>>& expr) {
//   MultMatMatLapack(expr.A(), expr.B(), C);
//   return A;
//
}  // namespace Tombino_bla
#endif
