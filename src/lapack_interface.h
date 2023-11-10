#ifndef FILE_LAPACK_INTERFACE_H
#define FILE_LAPACK_INTERFACE_H

#include <complex>
#include <iostream>
#include <string>

#include "matrix.h"
#include "vector.h"

// include std vector
#include <vector>

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

namespace ASC_bla
{

  // BLAS-1 functions:

  /*
    int daxpy_(integer *n, doublereal *da, doublereal *dx, 
    integer *incx, doublereal *dy, integer *incy);
  */
  // y += alpha x
  template <typename SX, typename SY>
  void AddVectorLapack (double alpha, VectorView<double,SX> x, VectorView<double,SY> y)
  {
    integer n = x.Size();
    integer incx = x.Dist();
    integer incy = y.Dist();
    int err = 
      daxpy_ (&n, &alpha, &x(0),  &incx, &y(0), &incy);
    if (err != 0)
      std::cerr << "daxpy returned errcode " << err << std::endl;
  }

  // BLAS-2 functions:

  // BLAS-3 functions:

  /*
    int dgemm_ (char *transa, char *transb, integer *m, integer *
    n, integer *k, doublereal *alpha, doublereal *a, integer *lda,
    doublereal *b, integer *ldb, doublereal *beta, doublereal *c__,
    integer *ldc);
  */

  // c = a*b
  template <ORDERING OA, ORDERING OB>
  void MultMatMatLapack(MatrixView<double, OA> a, MatrixView<double, OB> b,
                        MatrixView<double, ColMajor> c) {
    char transa_ = (OA == ColMajor) ? 'N' : 'T';
    char transb_ = (OB == ColMajor) ? 'N' : 'T';

    integer n = c.SizeRows();
    integer m = c.SizeCols();
    integer k = a.SizeCols();

    double alpha = 1.0;
    double beta = 0;
    integer lda = std::max(a.Dist(), 1ul);
    integer ldb = std::max(b.Dist(), 1ul);
    integer ldc = std::max(c.Dist(), 1ul);

    int err = dgemm_(&transa_, &transb_, &n, &m, &k, &alpha, a.Data(), &lda,
                     b.Data(), &ldb, &beta, c.Data(), &ldc);

    if (err != 0)
      throw std::runtime_error(std::string("MultMatMat got error") +
                               std::to_string(err));
    // MultMatMat got error"+std::to_string(err)));
  };

  template <ORDERING OA, ORDERING OB>
  void MultMatMatLapack(MatrixView<double, OA> a, MatrixView<double, OB> b,
                        MatrixView<double, RowMajor> c) {
    MultMatMatLapack(Transpose(b), Transpose(a), Transpose(c));
  };

  

template < ORDERING ORD>
  class LapackView
  {
    protected:
    Matrix<double, ORD> a;
   public:

  };

template <ORDERING ORD>
class LapackView : public LapackLU<LapackView<ORD>> {

 protected:
    Matrix<double, ORD> a;
    std::vector<integer> ipiv;

   public:
    LapackLU(Matrix<double, ORD> _a) : a(std::move(_a)), ipiv(a.SizeRows()) {
      integer m = a.SizeRows();
      if (m == 0) return;
      integer n = a.SizeCols();
      integer lda = a.Dist();
      integer info;

      // int dgetrf_(integer *m, integer *n, doublereal *a,
      //             integer * lda, integer *ipiv, integer *info);

      dgetrf_(&n, &m, &a(0, 0), &lda, &ipiv[0], &info);
    }

    // b overwritten with A^{-1} b
    template <typename Db>
    void Solve(VectorView<double, Db> b) const {
      char transa = (ORD == ColMajor) ? 'N' : 'T';
      integer n = a.SizeRows();
      integer nrhs = 1;
      integer lda = a.Dist();
      integer ldb = b.Size();
      integer info;

      // int dgetrs_(char *trans, integer *n, integer *nrhs,
      //             doublereal *a, integer *lda, integer *ipiv,
      //             doublereal *b, integer *ldb, integer *info);

      dgetrs_(&transa, &n, &nrhs, a.Data(), &lda, (integer*)ipiv.data(),
              b.Data(), &ldb, &info);
    }

    Matrix<double, ORD> Inverse() && {
      double hwork;
      integer lwork = -1;
      integer n = a.SizeRows();
      integer lda = a.Dist();
      integer info;

      // int dgetri_(integer *n, doublereal *a, integer *lda,
      //             integer *ipiv, doublereal *work, integer *lwork,
      //             integer *info);

      // query work-size
      dgetri_(&n, &a(0, 0), &lda, ipiv.data(), &hwork, &lwork, &info);
      lwork = integer(hwork);
      std::vector<double> work(lwork);
      dgetri_(&n, &a(0, 0), &lda, ipiv.data(), &work[0], &lwork, &info);
      return std::move(a);
    }

template <typename T, ORDERING ORD>
std::ostream& operator<<(std::ostream& os, const Matrix<T, ORD>& m) {
  os << std::endl;
  for (size_t i = 0; i < m.SizeRows(); i++) {
    for (size_t j = 0; j < m.SizeCols(); j++) {
      os << m(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}

    Matrix<double, ORD> LFactor() const {
      Matrix<double, ORD> L(a.SizeRows(), a.SizeCols());
      for (size_t i = 0; i < a.SizeRows(); i++) {
        for (size_t j = 0; j < a.SizeCols(); j++) {
          if (i > j)
            L(i, j) = a(i, j);
          else if (i == j)
            L(i, j) = 1;
          else
            L(i, j) = 0;
        }
      }
      return L;
    }

    Matrix<double, ORD> UFactor() const {
      Matrix<double, ORD> U(a.SizeRows(), a.SizeCols());
      for (size_t i = 0; i < a.SizeRows(); i++) {
        for (size_t j = 0; j < a.SizeCols(); j++) {
          if (i <= j)
            U(i, j) = a(i, j);
          else
            U(i, j) = 0;
        }
      }
      return U;
    }

    Matrix<double, ORD> PFactor() const {
      Matrix<double, ORD> P(a.SizeRows(), a.SizeCols());
      for (size_t i = 0; i < a.SizeRows(); i++) {
        for (size_t j = 0; j < a.SizeCols(); j++) {
          if (i == ipiv[j])
            P(i, j) = 1;
          else
            P(i, j) = 0;
        }
      }
      return P;
    };
  };

  template <ORDERING ORD>
  class LapackQR {
    //  QR decomposition using Hausholder reflection.
    //  A = Q R
    //  Q is orthogonal, R is upper triangular
    //  and Q is stored in the lower part of A as product of Hausholder
    //  reflections. R is stored in the upper part of A.

    Matrix<double, ORD> a;
    std::vector<double> tau;

   public:
    LapackQR(Matrix<double, ORD> _a) : a(std::move(_a)), tau(a.SizeCols()) {
      integer m = a.SizeRows();
      if (m == 0) return;
      integer n = a.SizeCols();
      integer lda = a.Dist();
      integer info;

      // int dgeqrf_(integer *m, integer *n, doublereal *a,
      //             integer *lda, doublereal *tau, doublereal *work,
      //             integer *lwork, integer *info);

      dgeqrf_(&n, &m, &a(0, 0), &lda, &tau[0], &tau[0], &info);
      // the upper triangular part of a contains R
      // the lower triangular part of a contains the Householder vectors.

      if (info != 0)
        throw std::runtime_error(std::string("LapackQR got error") +
                                 std::to_string(info));
    }

    // b overwritten with A^{-1} b
    template <typename Db>
    void Solve(VectorView<double, Db> b) const {
      char transa = (ORD == ColMajor) ? 'N' : 'T';
      integer n = a.SizeRows();
      integer nrhs = 1;
      integer lda = a.Dist();
      integer ldb = b.Size();
      integer info;

      // int dgetrs_(char *trans, integer *n, integer *nrhs,
      //             doublereal *a, integer *lda, integer *ipiv,
      //             doublereal *b, integer *ldb, integer *info);

      dgetrs_(&transa, &n, &nrhs, a.Data(), &lda, (integer*)tau.data(),
              b.Data(), &ldb, &info);
    };

    Matrix<double, ORD> Inverse() && {
      double hwork;
      integer lwork = -1;
      integer n = a.SizeRows();
      integer lda = a.Dist();
      integer info;

      // int dorgqr_(integer *m, integer *n, integer *k,
      //             doublereal *a, integer *lda, doublereal *tau,
      //             doublereal *work, integer *lwork, integer *info);

      // query work-size
      dorgqr_(&n, &n, &n, &a(0, 0), &lda, &tau[0], &hwork, &lwork, &info);
      lwork = integer(hwork);
      std::vector<double> work(lwork);
      dorgqr_(&n, &n, &n, &a(0, 0), &lda, &tau[0], &work[0], &lwork, &info);
      return std::move(a);
    };

    Matrix<double, ORD> RFactor() const {
      Matrix<double, ORD> R(a.SizeRows(), a.SizeCols());
      for (size_t i = 0; i < a.SizeRows(); i++) {
        for (size_t j = 0; j < a.SizeCols(); j++) {
          if (i <= j)
            R(i, j) = a(i, j);
          else
            R(i, j) = 0;
        }
      }
      return R;
    };

    // to get the Q factor we need to apply the Householder reflections with the
    // function dorgqr
    Matrix<double, ORD> QFactor() const {
      Matrix<double, ORD> Q(a.SizeRows(), a.SizeCols());
      for (size_t i = 0; i < a.SizeRows(); i++) {
        for (size_t j = 0; j < a.SizeCols(); j++) {
          if (i == j)
            Q(i, j) = 1;
          else
            Q(i, j) = 0;
        }
      }
      integer m = a.SizeRows();
      integer n = a.SizeCols();
      integer lda = a.Dist();
      integer lwork = -1;
      double hwork;
      integer info;
      dorgqr_(&n, &m, &n, &Q(0, 0), &lda, &tau[0], &hwork, &lwork, &info);
      lwork = integer(hwork);
      std::vector<double> work(lwork);
      dorgqr_(&n, &m, &n, &Q(0, 0), &lda, &tau[0], &work[0], &lwork, &info);
      return Q;
    };
  };

  // create the class LapackEVP to calculate the eigenvalues and eigenvectors
  // using dgeev to calculate the right eignevalues and eigenvectors
  template <ORDERING ORD>
  class LapackEVP {
    Matrix<double, ORD> a;
    std::vector<double> wr;
    std::vector<double> wi;
    std::vector<double> vl;
    std::vector<double> vr;

   public:
    LapackEVP(Matrix<double, ORD> _a)
        : a(std::move(_a)),
          wr(a.SizeCols()),
          wi(a.SizeCols()),
          vl(a.SizeCols() * a.SizeCols()),
          vr(a.SizeCols() * a.SizeCols()) {
      integer m = a.SizeRows();
      if (m == 0) return;
      integer n = a.SizeCols();
      integer lda = a.Dist();
      integer info;

      // int dgeev_(char *jobvl, char *jobvr, integer *n, doublereal *a,
      //            integer *lda, doublereal *wr, doublereal *wi,
      //            doublereal *vl, integer *ldvl, doublereal *vr,
      //            integer *ldvr, doublereal *work, integer *lwork,
      //            integer *info);

      dgeev_("N", "V", &n, &m, &a(0, 0), &lda, &wr[0], &wi[0], &vl[0], &n,
             &vr[0], &n, &vl[0], &n, &info);
      // the upper triangular part of a contains R
      // the lower triangular part of a contains the Householder vectors.

      if (info != 0)
        throw std::runtime_error(std::string("LapackEVP got error") +
                                 std::to_string(info));
    }

    // b overwritten with A^{-1} b
  };

  // create the class LapackSVD to calculate the eigenvalues and eigenvectors
  // using the SVD method
  template <ORDERING ORD>
  class LapackSVD {
    Matrix<double, ORD> a;
    std::vector<double> s;
    std::vector<double> u;
    std::vector<double> vt;

   public:
    LapackSVD(Matrix<double, ORD> _a)
        : a(std::move(_a)),
          s(std::min(a.SizeCols(), a.SizeRows())),
          u(a.SizeCols() * a.SizeCols()),
          vt(a.SizeRows() * a.SizeRows()) {
      integer m = a.SizeRows();
      if (m == 0) return;
      integer n = a.SizeCols();
      integer lda = a.Dist();
      integer info;

      // int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
      //             doublereal *a, integer *lda, doublereal *s,
      //             doublereal *u, integer *ldu, doublereal *vt,
      //             integer *ldvt, doublereal *work, integer *lwork,
      //             integer *info);

      dgesvd_("A", "A", &n, &m, &a(0, 0), &lda, &s[0], &u[0], &n, &vt[0], &n,
              &u[0], &n, &info);
      // the upper triangular part of a contains R
      // the lower triangular part of a contains the Householder vectors.

      if (info != 0)
        throw std::runtime_error(std::string("LapackSVD got error") +
                                 std::to_string(info));
    }

    // b overwritten with A^{-1} b
  };

// lapack view
template < ORDERING ORD>
std::ostream& operator<<(std::ostream& os, const LapackView< ORD>& lpv) {
  os << std::endl;
  Matrix<double, ORD> m (lpv.Mat())
  for (size_t i = 0; i < lpv.SizeRows(); i++) {
    for (size_t j = 0; j < m.SizeCols(); j++) {
      os << m(i, j) << " ";
    }
    os << std::endl;
  }
  return os;
}
  }  // namespace ASC_bla

#endif
