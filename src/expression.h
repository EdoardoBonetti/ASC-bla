#ifndef FILE_EXPRESSION_H
#define FILE_EXPRESSION_H
/**************************************************************************/
/* File:   matrix.h                                                       */
/* Author: Joachim Schoeberl, Edoardo Bonetti                             */
/* Date:   13. 10. 2023                                                   */
/**************************************************************************/
#include <complex>
#include <concepts>
#include <iostream>
#include <type_traits>

// #include "lapack_interface.h"

namespace Tombino_bla {

typedef std::complex<double> dcomplex;

template <typename T>
concept ValidSCAL = std::is_same_v<T, int> || std::is_same_v<T, double> || std::is_same_v<T, dcomplex>;

/****************************************************************/
/*                   Expressions for Vectors                    */
/****************************************************************/

// create a template function that can be used only with double or complex

template <typename T>
class VecExpr {
 public:
  auto Upcast() const { return static_cast<const T&>(*this); }
  size_t Size() const { return Upcast().Size(); }
  size_t Dist() const { return Upcast().Dist(); }
  auto operator()(size_t i) const { return Upcast()(i); }
  auto& operator()(size_t i) { return Upcast()(i); }
};

template <typename TA, typename TB>
class SumVecExpr : public VecExpr<SumVecExpr<TA, TB>> {
  TA a_;
  TB b_;

 public:
  SumVecExpr(TA a, TB b) : a_(a), b_(b) {
    if (a_.Size() != b_.Size()) {
      std::cout << "Error: the two vectors have different sizes" << std::endl;
      exit(1);
    }
  }
  auto operator()(size_t i) const { return a_(i) + b_(i); }
  size_t Size() const { return a_.Size(); }
  size_t Dist() const { return a_.Dist(); }
};

template <typename TA, typename TB>
auto operator+(const VecExpr<TA>& a, const VecExpr<TB>& b) {
  return SumVecExpr(a.Upcast(), b.Upcast());
}

template <typename TA, typename TB>
auto operator+=(const VecExpr<TA>& a, const VecExpr<TB>& b)
{
  return SumVecExpr(a.Upcast(), b.Upcast());
}

template <typename TA, typename TB>
class SubVecExpr : public VecExpr<SubVecExpr<TA, TB>> {
  TA a_;
  TB b_;

 public:
  SubVecExpr(TA a, TB b) : a_(a), b_(b) {
    if (a_.Size() != b_.Size()) {
      std::cout << "Error: the two vectors have different sizes" << std::endl;
      exit(1);
    }
  }
  auto operator()(size_t i) const { return a_(i) - b_(i); }
  size_t Size() const { return a_.Size(); }
  size_t Dist() const { return a_.Dist(); }
};

template <typename TA, typename TB>
auto operator-(const VecExpr<TA>& a, const VecExpr<TB>& b) {
  return SubVecExpr(a.Upcast(), b.Upcast());
}

template <typename TA, typename TB>
auto operator-=(const VecExpr<TA>& a, const VecExpr<TB>& b)
{
  return SubVecExpr(a.Upcast(), b.Upcast());
}

template <typename TSCAL, typename TV>
class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL, TV>> {
  TSCAL scal_;
  TV vec_;

 public:
  ScaleVecExpr(TSCAL scal, TV vec) : scal_(scal), vec_(vec) {}
  auto operator()(size_t i) const { return scal_ * vec_(i); }
  size_t Size() const { return vec_.Size(); }
  size_t Dist() const { return vec_.Dist(); }
};
template <ValidSCAL TSCAL, typename T>
auto operator*(TSCAL scal, const VecExpr<T>& v) {
  return ScaleVecExpr(scal, v.Upcast());
}

template <ValidSCAL TSCAL, typename T>
auto operator*=(const VecExpr<T>& v, TSCAL scal)
{
  return ScaleVecExpr(1 + scal, v.Upcast());
}

template <typename TSCAL, typename TV>
class ScaleVecSumExpr : public VecExpr<ScaleVecSumExpr<TSCAL, TV>> {
  TSCAL scal_;
  TV vec_;

 public:
  ScaleVecSumExpr(TSCAL scal, TV vec) : scal_(scal), vec_(vec) {}
  auto operator()(size_t i) const { return scal_ + vec_(i); }
  size_t Size() const { return vec_.Size(); }
  size_t Dist() const { return vec_.Dist(); }
};

template <ValidSCAL TSCAL, typename T>
auto operator+(TSCAL scal, const VecExpr<T>& v) {
  return ScaleVecSumExpr(scal, v.Upcast());
}

template <typename TSCAL, typename TV>
class ScaleVecSubExpr : public VecExpr<ScaleVecSubExpr<TSCAL, TV>> {
  TSCAL scal_;
  TV vec_;

 public:
  ScaleVecSubExpr(TSCAL scal, TV vec) : scal_(scal), vec_(vec) {}
  auto operator()(size_t i) const { return scal_ - vec_(i); }
  size_t Size() const { return vec_.Size(); }
  size_t Dist() const { return vec_.Dist(); }
};

template <ValidSCAL TSCAL, typename T>
auto operator-(TSCAL scal, const VecExpr<T>& v) {
  return ScaleVecSubExpr(scal, v.Upcast());
}

template <typename TA, typename TB>
auto InnerProduct(const VecExpr<TA>& a, const VecExpr<TB>& b) {
  typedef decltype(std::declval<TA>()(0) * std::declval<TB>()(0)) TRES;
  TRES sum = a(0) * b(0);
  for (size_t i = 0; i < a.Size(); i++) sum = sum + (a(i) * b(i));
  return sum;
}

template <typename TA, typename TB>
auto operator*(const VecExpr<TA>& a, const VecExpr<TB>& b) {
  if (a.Size() != b.Size()) {
    std::cout << "Error: the two vectors have different sizes" << std::endl;
    exit(1);
  }
  return InnerProduct(a.Upcast(), b.Upcast());
}

// L2 norm of a vector
template <typename T>
auto L2Norm(const VecExpr<T>& v)
{
  return std::sqrt(InnerProduct(v, v));
}

/****************************************************************/
/*                   Expressions for Matrices                    */
/****************************************************************/

template <typename T>
class MatExpr {
 public:
  auto Upcast() const { return static_cast<const T&>(*this); }
  size_t SizeCols() const { return Upcast().SizeCols(); }
  size_t SizeRows() const { return Upcast().SizeRows(); }
  auto operator()(size_t i, size_t j) const { return Upcast()(i, j); }
  auto& operator()(size_t i, size_t j) { return Upcast()(i, j); }
};

template <typename TA, typename TB>
class SumMatExpr : public MatExpr<SumMatExpr<TA, TB>> {
  TA a_;
  TB b_;

 public:
  SumMatExpr(TA a, TB b) : a_(a), b_(b) {
    if ((a_.SizeCols() != b_.SizeCols()) || (a_.SizeRows() != b_.SizeRows())) {
      std::cout << "Error: the two matrices have different sizes" << std::endl;
      exit(1);
    }
  }
  auto operator()(size_t i, size_t j) const { return a_(i, j) + b_(i, j); }
  size_t SizeCols() const { return a_.SizeCols(); }
  size_t SizeRows() const { return a_.SizeRows(); }
};

template <typename TA, typename TB>
auto operator+(const MatExpr<TA>& a, const MatExpr<TB>& b) {
  return SumMatExpr(a.Upcast(), b.Upcast());
}

// operator += for matrices
template <typename TA, typename TB>
auto operator+=(const MatExpr<TA>& a, const MatExpr<TB>& b)
{
  return SumMatExpr(a.Upcast(), b.Upcast());
}

template <typename TA, typename TB>
class SubMatExpr : public MatExpr<SubMatExpr<TA, TB>> {
  TA a_;
  TB b_;

 public:
  SubMatExpr(TA a, TB b) : a_(a), b_(b) {
    if ((a_.SizeCols() != b_.SizeCols()) || (a_.SizeRows() != b_.SizeRows())) {
      std::cout << "Error: the two matrices have different sizes" << std::endl;
      exit(1);
    }
  }
  auto operator()(size_t i, size_t j) const { return a_(i, j) - b_(i, j); }
  size_t SizeCols() const { return a_.SizeCols(); }
  size_t SizeRows() const { return a_.SizeRows(); }
};

template <typename TA, typename TB>
auto operator-(const MatExpr<TA>& a, const MatExpr<TB>& b) {
  return SubMatExpr(a.Upcast(), b.Upcast());
}
template <typename TA, typename TB>
auto operator-=(const MatExpr<TA>& a, const MatExpr<TB>& b)
{
  return SubMatExpr(a.Upcast(), b.Upcast());
}

template <typename TA, typename TB>
class ProdMatExpr : public MatExpr<ProdMatExpr<TA, TB>> {
  TA a_;
  TB b_;

 public:
  ProdMatExpr(TA a, TB b) : a_(a), b_(b) {
    if (a_.SizeCols() != b_.SizeRows()) {
      std::cout << "Error: the two matrices cannot be multiplied" << std::endl;
      exit(1);
    }
  }
  auto operator()(size_t i, size_t j) const {
    typedef decltype(std::declval<TA>()(0, 0) * std::declval<TB>()(0, 0)) TRES;
    TRES sum = a_(i, 0) * b_(0, j);
    for (size_t k = 1; k < a_.SizeCols(); k++) {
      sum += a_(i, k) * b_(k, j);
    }
    return sum;
  }
  size_t SizeCols() const { return b_.SizeCols(); }
  size_t SizeRows() const { return a_.SizeRows(); }
};

template <typename TA, typename TB>
auto operator*(const MatExpr<TA>& a, const MatExpr<TB>& b) {
  return ProdMatExpr(a.Upcast(), b.Upcast());
}

template <typename TSCAL, typename TM>
class ScaleMatExpr : public MatExpr<ScaleMatExpr<TSCAL, TM>>
{
  TSCAL scal_;
  TM mat_;

 public:
  ScaleMatExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
  auto Upcast() { return ScalMatExpr(scal_, mat_); }
  auto operator()(size_t i, size_t j) const { return scal_ * mat_(i, j); }
  size_t SizeCols() const { return mat_.SizeCols(); }
  size_t SizeRows() const { return mat_.SizeRows(); }
};

template <ValidSCAL TSCAL, typename T>
auto operator*(TSCAL scal, const MatExpr<T>& m)
{
  return ScaleMatExpr(scal, m.Upcast());
}

template <ValidSCAL TSCAL, typename T>
auto operator*=(const MatExpr<T>& m, TSCAL scal)
{
  return ScaleMatExpr(1 + scal, m.Upcast());
}

template <typename TSCAL, typename TM>
class ScaleMatSumExpr : public MatExpr<ScaleMatSumExpr<TSCAL, TM>>
{
  TSCAL scal_;
  TM mat_;

 public:
  ScaleMatSumExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
  auto Upcast() { return ScaleMatSumExpr(scal_, mat_); }
  auto operator()(size_t i, size_t j) const { return scal_ + mat_(i, j); }
  size_t SizeCols() const { return mat_.SizeCols(); }
  size_t SizeRows() const { return mat_.SizeRows(); }
};

template <ValidSCAL TSCAL, typename T>
auto operator+(TSCAL scal, const MatExpr<T>& m)
{
  return ScaleMatSumExpr(scal, m.Upcast());
}

template <ValidSCAL TSCAL, typename T>
auto operator+=(const MatExpr<T>& m, TSCAL scal)
{
  return ScaleMatSumExpr(1 + scal, m.Upcast());
}
// define a matrix scalar subtraction
template <typename TSCAL, typename TM>
class ScaleMatSubExpr : public MatExpr<ScaleMatSubExpr<TSCAL, TM>>
{
  TSCAL scal_;
  TM mat_;

 public:
  ScaleMatSubExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
  auto Upcast() { return ScaleMatSubExpr(scal_, mat_); }
  auto operator()(size_t i, size_t j) const { return scal_ - mat_(i, j); }
  size_t SizeCols() const { return mat_.SizeCols(); }
  size_t SizeRows() const { return mat_.SizeRows(); }
};

template <ValidSCAL TSCAL, typename T>
auto operator-(TSCAL scal, const MatExpr<T>& m)
{
  return ScaleMatSubExpr(scal, m.Upcast());
}

template <ValidSCAL TSCAL, typename T>
auto operator-=(const MatExpr<T>& m, TSCAL scal)
{
  return ScaleMatSubExpr(1 + scal, m.Upcast());
}

/****************************************************************/
/*                   Mixed Expressions Vectors                    */
/****************************************************************/

template <typename TM, typename TV>
class MatVecExpr : public VecExpr<MatVecExpr<TM, TV>>
{
  TM mat_;
  TV vec_;

 public:
  MatVecExpr(TM mat, TV vec) : mat_(mat), vec_(vec)
  {
    if (mat_.SizeCols() != vec_.Size())
    {
      std::cout << "Error: the matrix and the vector cannot be multiplied"
                << std::endl;
      exit(1);
    }
  }
  auto Upcast() { return MatVecExpr(mat_, vec_); }
  auto operator()(size_t i) const
  {
    typedef decltype(std::declval<TM>()(0, 0) * std::declval<TV>()(0)) TRES;
    TRES sum = mat_(i, 0) * vec_(0);
    for (size_t j = 1; j < mat_.SizeCols(); j++)
    {
      sum += mat_(i, j) * vec_(j);
    }
    return sum;
  }
  size_t Size() const { return mat_.SizeRows(); }
  size_t Dist() const { return vec_.Dist(); }
};

template <typename TM, typename TV>
auto operator*(const MatExpr<TM>& m, const VecExpr<TV>& v)
{
  return MatVecExpr(m.Upcast(), v.Upcast());
}

template <typename TM, typename TV>
auto operator*=(const VecExpr<TV>& v, const MatExpr<TM>& m)
{
  return MatVecExpr(1 + m.Upcast(), v.Upcast());
}

template <typename VA, typename VB>
class OuterVecVec : public MatExpr<OuterVecVec<VA, VB>> {
  VA a_;
  VB b_;

 public:
  OuterVecVec(VA a, VB b) : a_(a), b_(b) {}
  auto Upcast() { return OuterVecVec(a_, b_); }
  auto operator()(size_t i, size_t j) const {
    typedef decltype(std::declval<VA>()(i) * std::declval<VB>()(j)) TRES;
    TRES sum = a_(i) * b_(j);
    return sum;
  }
  size_t SizeCols() const { return b_.Size(); }
  size_t SizeRows() const { return a_.Size(); }
};

template <typename VA, typename VB>
auto operator^(const VecExpr<VA>& a, const VecExpr<VB>& b) {
  return OuterVecVec(a.Upcast(), b.Upcast());
}

template <typename T>
std::ostream& operator<<(std::ostream& ost, const VecExpr<T>& v) {
  if (v.Size() > 0) ost << v(0);
  for (size_t i = 1; i < v.Size(); i++) ost << ", " << v(i);
  return ost;
}

template <typename T>
std::ostream& operator<<(std::ostream& ost, const MatExpr<T>& m) {
  for (size_t i = 0; i < m.SizeRows(); i++) {
    for (size_t j = 0; j < m.SizeCols(); j++) {
      ost << m(i, j) << " ";
    }
    ost << std::endl;
  }
  return ost;
}

/****************************************************************/
/*                   Expressions for Lapack                     */
/****************************************************************/

// now we create a syntax for the lapack interface

/****************************************************************/
/*                   Expressions for Tensors                    */
/****************************************************************/
// need to implement the tensor expression, in the tensor expression we need to
// define the following operations:
// 1) tensor + tensor
// 2) tensortensor contraction

}  // namespace Tombino_bla

#endif
