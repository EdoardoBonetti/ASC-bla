#ifndef FILE_EXPRESSION_H
#define FILE_EXPRESSION_H
/**************************************************************************/
/* File:   matrix.h                                                       */
/* Author: Joachim Schoeberl, Edoardo Bonetti                             */
/* Date:   13. 10. 2023                                                   */
/**************************************************************************/
namespace ASC_bla
{

  template <typename T>
  class VecExpr
  {
   public:
    auto Upcast() const { return static_cast<const T&>(*this); }
    size_t Size() const { return Upcast().Size(); }
    size_t Dist() const { return Upcast().Dist(); }
    auto operator()(size_t i) const { return Upcast()(i); }
    auto& operator()(size_t i) { return Upcast()(i); }
  };
  
 
  template <typename TA, typename TB>
  class SumVecExpr : public VecExpr<SumVecExpr<TA,TB>>
  {
    TA a_;
    TB b_;
  public:
    SumVecExpr (TA a, TB b) : a_(a), b_(b) { }
    auto Upcast() { return SumExpr(a_, b_); }
    auto operator() (size_t i) const { return a_(i)+b_(i); }
    size_t Size() const { return a_.Size(); }
    size_t Dist() const { return a_.Dist(); }
  };
  
  template <typename TA, typename TB>
  auto operator+ (const VecExpr<TA> & a, const VecExpr<TB> & b)
  {
    return SumVecExpr(a.Upcast(), b.Upcast());
  }

  template <typename TSCAL, typename TV>
  class ScaleVecExpr : public VecExpr<ScaleVecExpr<TSCAL, TV>> {
    TSCAL scal_;
    TV vec_;
  public:
    ScaleVecExpr (TSCAL scal, TV vec) : scal_(scal), vec_(vec) { }
    auto Upcast() { return ScalVecExpr(scal_, vec_); }
    auto operator() (size_t i) const { return scal_*vec_(i); }
    size_t Size() const { return vec_.Size(); }
    size_t Dist() const { return vec_.Dist(); }
  };

  template <typename T>
  auto operator* (double scal, const VecExpr<T> & v)
  {
    return ScaleVecExpr(scal, v.Upcast());
  }

  // addition by a scalar

  template <typename TSCAL, typename TV>
  class ScaleVecSumExpr : public VecExpr<ScaleVecSumExpr<TSCAL, TV>> {
    TSCAL scal_;
    TV vec_;

   public:
    ScaleVecSumExpr(TSCAL scal, TV vec) : scal_(scal), vec_(vec) {}
    auto Upcast() { return ScalVecExpr(scal_, vec_); }
    auto operator()(size_t i) const { return scal_ + vec_(i); }
    size_t Size() const { return vec_.Size(); }
    size_t Dist() const { return vec_.Dist(); }
  };

  template <typename T>
  auto operator+(double scal, const VecExpr<T>& v) {
    return ScaleVecSumExpr(scal, v.Upcast());
  }

  /****************************************************************/
  /*         My implementation of the matrix class                */
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
    SumMatExpr(TA a, TB b) : a_(a), b_(b) {}
    auto Upcast() { return SumExpr(a_, b_); }
    auto operator()(size_t i, size_t j) const { return a_(i, j) + b_(i, j); }
    size_t SizeCols() const { return a_.SizeCols(); }
    size_t SizeRows() const { return a_.SizeRows(); }
  };

  template <typename TA, typename TB>
  auto operator+(const MatExpr<TA>& a, const MatExpr<TB>& b) {
    return SumMatExpr(a.Upcast(), b.Upcast());
    // return SumMatExpr<TA, TB>{a.Upcast(), b.Upcast()};
  }

  template <typename TSCAL, typename TM>
  class ScaleMatExpr : public MatExpr<ScaleMatExpr<TSCAL, TM>> {
    TSCAL scal_;
    TM mat_;

   public:
    ScaleMatExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
    auto Upcast() { return ScalMatExpr(scal_, mat_); }
    auto operator()(size_t i, size_t j) const { return scal_ * mat_(i, j); }
    size_t SizeCols() const { return mat_.SizeCols(); }
    size_t SizeRows() const { return mat_.SizeRows(); }
  };

  template <typename T>
  auto operator*(double scal, const MatExpr<T>& m) {
    // return ScaleMatExpr<double, T>(scal, m.Upcast());
    //  Why not ?
    return ScaleMatExpr(scal, m.Upcast());
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

  // define a matrix vector product
  template <typename TM, typename TV>
  class MatVecExpr : public VecExpr<MatVecExpr<TM, TV>> {
    TM mat_;
    TV vec_;

   public:
    MatVecExpr(TM mat, TV vec) : mat_(mat), vec_(vec) {}
    auto Upcast() { return MatVecExpr(mat_, vec_); }
    auto operator()(size_t i) const {
      auto sum = 0.0;
      for (size_t j = 0; j < mat_.SizeCols(); j++) {
        sum += mat_(i, j) * vec_(j);
      }
      return sum;
    }
    size_t Size() const { return mat_.SizeRows(); }
    size_t Dist() const { return vec_.Dist(); }
  };

  template <typename TM, typename TV>
  auto operator*(const MatExpr<TM>& m, const VecExpr<TV>& v) {
    return MatVecExpr(m.Upcast(), v.Upcast());
  }

  // define a matrix vector product
  template <typename TA, typename TB>
  class ProdMatExpr : public MatExpr<ProdMatExpr<TA, TB>> {
    TA a_;
    TB b_;

   public:
    ProdMatExpr(TA a, TB b) : a_(a), b_(b) {}
    auto Upcast() { return ProdMatExpr(a_, b_); }
    auto operator()(size_t i, size_t j) const {
      auto sum = 0.0;
      for (size_t k = 0; k < a_.SizeCols(); k++) {
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

  // define a matrix scalar addition
  template <typename TSCAL, typename TM>
  class ScaleMatSumExpr : public MatExpr<ScaleMatSumExpr<TSCAL, TM>> {
    TSCAL scal_;
    TM mat_;

   public:
    ScaleMatSumExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
    auto Upcast() { return ScaleMatSumExpr(scal_, mat_); }
    auto operator()(size_t i, size_t j) const { return scal_ + mat_(i, j); }
    size_t SizeCols() const { return mat_.SizeCols(); }
    size_t SizeRows() const { return mat_.SizeRows(); }
  };

  template <typename T>
  auto operator+(double scal, const MatExpr<T>& m) {
    return ScaleMatSumExpr(scal, m.Upcast());
  }

  // define a matrix scalar subtraction
  template <typename TSCAL, typename TM>
  class ScaleMatSubExpr : public MatExpr<ScaleMatSubExpr<TSCAL, TM>> {
    TSCAL scal_;
    TM mat_;

   public:
    ScaleMatSubExpr(TSCAL scal, TM mat) : scal_(scal), mat_(mat) {}
    auto Upcast() { return ScaleMatSubExpr(scal_, mat_); }
    auto operator()(size_t i, size_t j) const { return scal_ - mat_(i, j); }
    size_t SizeCols() const { return mat_.SizeCols(); }
    size_t SizeRows() const { return mat_.SizeRows(); }
  };

  template <typename T>
  auto operator-(double scal, const MatExpr<T>& m) {
    return ScaleMatSubExpr(scal, m.Upcast());
  }

  // define a matrix scalar subtraction
  template <typename TA, typename TB>
  class SumMatSubExpr : public MatExpr<SumMatSubExpr<TA, TB>> {
    TA a_;
    TB b_;

   public:
    SumMatSubExpr(TA a, TB b) : a_(a), b_(b) {}
    auto Upcast() { return SumMatSubExpr(a_, b_); }
    auto operator()(size_t i, size_t j) const { return a_(i, j) - b_(i, j); }
    size_t SizeCols() const { return a_.SizeCols(); }
    size_t SizeRows() const { return a_.SizeRows(); }
  };

  template <typename TA, typename TB>
  auto operator-(const MatExpr<TA>& a, const MatExpr<TB>& b) {
    return SumMatSubExpr(a.Upcast(), b.Upcast());
  }

  // define and inner product between two vectors
  template <typename TA, typename TB>
  auto operator*(const VecExpr<TA>& a, const VecExpr<TB>& b) {
    return InnerProduct(a.Upcast(), b.Upcast());
  }

  template <typename TA, typename TB>
  auto InnerProduct(const VecExpr<TA>& a, const VecExpr<TB>& b) {
    auto sum = 0.0;
    for (size_t i = 0; i < a.Size(); i++) {
      sum += a(i) * b(i);
    }
    return sum;
  }
}
 
#endif
