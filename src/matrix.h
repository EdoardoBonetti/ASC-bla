#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"
#include "vector.h"

namespace Tombino_bla {
enum ORDERING { RowMajor, ColMajor };

template <typename T = double, ORDERING ORD = ORDERING::RowMajor>
class MatrixView : public MatExpr<MatrixView<T, ORD>> {
 protected:
  T* data_;
  size_t rows_;
  size_t cols_;
  size_t d_r_;
  size_t d_c_;

 public:
  // Constructor
  MatrixView(size_t rows, size_t cols, T* data)
      : data_(data), rows_(rows), cols_(cols) {
    if constexpr (ORD == RowMajor) {
      d_c_ = 1;
      d_r_ = cols;
    } else {
      d_c_ = rows;
      d_r_ = 1;
    }
  }

  MatrixView(size_t rows, size_t cols, size_t d_r, size_t d_c, T* data)
      : data_(data), rows_(rows), cols_(cols), d_r_(d_r), d_c_(d_c) {}

  // Copy constructor from MatrixView
  MatrixView& operator=(const MatrixView& m2) {
    *this = static_cast<const MatExpr<MatrixView<T, ORD>>&>(m2);
    return *this;
  }
  // Copy assignment operator, needs to be modified for row maj
  template <typename TB>
  MatrixView& operator=(const MatExpr<TB>& m2) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) = m2(i, j);
      }
    }
    return *this;
  }

  MatrixView& operator=(T scal) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) = scal;
      }
    }
    return *this;
  }

  auto View() const { return MatrixView(rows_, cols_, d_r_, d_c_, data_); }
  size_t SizeCols() const { return cols_; }
  size_t SizeRows() const { return rows_; }
  T* Data() const { return data_; }
  size_t DistRows() const { return d_r_; }
  size_t DistCols() const { return d_c_; }
  T& operator()(size_t i, size_t j) { return data_[i * d_r_ + j * d_c_]; }
  const T& operator()(size_t i, size_t j) const {
    return data_[i * d_r_ + j * d_c_];
  }
  // T& operator()(size_t i) { return this->operator()(i); }
  // const T& operator()(size_t i) const { return this->operator()(i); }

  auto Row(size_t i) const {
    if constexpr (ORD == RowMajor) {
      // extract the i-th row: in particular the i-th row is a vector of size
      // cols_

      return VectorView<T>(cols_, data_ + i * d_r_);
    } else {
      return VectorView<T, size_t>(cols_, d_c_, data_ + i * d_r_);
    }
  }

  auto Col(size_t i) const {
    if constexpr (ORD == RowMajor) {
      return VectorView<T, size_t>(rows_, d_r_, data_ + i * d_c_);
    } else {
      return VectorView<T>(rows_, data_ + i * d_c_);
    }
  }

  auto Rows(size_t first, size_t next) const {
    // if constexpr (ORD == RowMajor) {
    //   return MatrixView<T, ORD>(next - first, cols_, d_r_, d_c_, data_ +
    //   first * d_r_);
    // } else {
    return MatrixView<T, ORD>(next - first, cols_, d_r_, d_c_,
                              data_ + first * d_r_);
    //}
  }

  auto Cols(size_t first, size_t next) const {
    // if constexpr (ORD == RowMajor) {
    //   return MatrixView<T, ORD>(rows_, next - first, d_r_, d_c_, data_ +
    //   first * d_c_);
    // } else {
    return MatrixView<T, ORD>(rows_, next - first, d_r_, d_c_,
                              data_ + first * d_c_);
    // }
  }

  // Ranges

  // Slices

  /*
  auto RowSlice(size_t i, size_t first, size_t slice) const {
    if constexpr (ORD == RowMajor) {
      std::cout << "i = " << i << " first = " << first << " slice = " << slice
  << std::endl;

      std::cout << "VectorView<T, size_t>(cols_ / slice, d_c_ * slice, data_ + i
  * d_r_ + first * d_c_);" << std::endl; std::cout << "numbs Elements " << cols_
  / slice << std::endl; std::cout << "jump " << d_c_ * slice << std::endl;
      std::cout << "position first element " << i * d_r_ + first * d_c_ <<
  std::endl;

      return VectorView<T, size_t>(cols_ / slice, d_c_ * slice, data_ + i * d_r_
  + first * d_c_); } else { return VectorView<T, size_t>(cols_ / slice, d_c_ *
  slice, data_ + i * d_r_ + first * d_c_);
    }
  };

  auto ColSlice(size_t j, size_t first, size_t slice) const {
    if constexpr (ORD == RowMajor) {
      return VectorView<T, size_t>(rows_ / slice, d_r_ * slice, data_ + j * d_c_
  + first * d_r_); } else { return VectorView<T, size_t>(rows_ / slice, d_r_ *
  slice, data_ + j * d_c_ + first * d_r_);
    }
  };
  */
  auto RSlice(size_t first, size_t slice) const {
    if constexpr (ORD == RowMajor) {
      return MatrixView<T, ORD>(rows_, (cols_ + 1) / slice, d_r_, d_c_ * slice,
                                data_ + first * d_c_);
    } else {
      return MatrixView<T, ORD>(rows_, cols_ / slice, d_r_, d_c_ * slice,
                                data_ + first * d_c_);
    }
  };

  auto CSlice(size_t first, size_t slice) const {
    if constexpr (ORD == RowMajor) {
      std::cout << "rows_ = " << rows_ << " cols_ = " << cols_ << std::endl;
      std::cout << "d_r_ = " << d_r_ << " d_c_ = " << d_c_ << std::endl;

      std::cout << "rows_ / slice = " << rows_ / slice << std::endl;

      return MatrixView<T, ORD>((rows_ + 1) / slice, cols_, d_r_ * slice, d_c_,
                                data_ + first * d_r_);
    } else {
      return MatrixView<T, ORD>(rows_ / slice, cols_, d_r_ * slice, d_c_,
                                data_ + first * d_r_);
    }
  };
  // now create a getter and setter for diagonal elements
  auto Diag() const {
    if constexpr (ORD == RowMajor) {
      return VectorView<T, size_t>(rows_, d_r_ + d_c_, data_);
    } else {
      return VectorView<T, size_t>(rows_, d_r_ + d_c_, data_);
    }
  };

  // operator+= mat
  template <typename TB>
  MatrixView& operator+=(const MatExpr<TB>& m2) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) += m2(i, j);
      }
    }
    return *this;
  }

  // operator-= mat
  template <typename TB>
  MatrixView& operator-=(const MatExpr<TB>& m2) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) -= m2(i, j);
      }
    }
    return *this;
  }

  // operator*= scal
  template <ValidSCAL TSCAL>
  MatrixView& operator*=(TSCAL scal) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) *= scal;
      }
    }
    return *this;
  }

  // operator/= scal
  template <ValidSCAL TSCAL>
  MatrixView& operator/=(TSCAL& scal) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        (*this)(i, j) /= scal;
      }
    }
    return *this;
  }
};

// template <typename T = double, ORDERING ORD = ORDERING::RowMajor>
// class MatrixView;

template <typename T = double, ORDERING ORD = ORDERING::RowMajor>
class Matrix;

template <typename T, ORDERING ORD>
class Matrix : public MatrixView<T, ORD> {
  typedef MatrixView<T, ORD> BASE;
  using BASE::cols_;
  using BASE::d_c_;
  using BASE::d_r_;
  using BASE::data_;
  using BASE::rows_;

 public:
  Matrix(size_t rows, size_t cols)
      : MatrixView<T, ORD>(rows, cols, new T[rows * cols]) {
    ;
  }
  Matrix(size_t size) : MatrixView<T, ORD>(size, size, new T[size * size]) { ; }
  // Copy constructor from Matrix
  Matrix(const Matrix& m) : Matrix(m.SizeRows(), m.SizeCols()) { *this = m; }

  // Copy constructor from MatrixView: New feature
  // Matrix(const MatrixView<T, ORD>& m) : Matrix(m.SizeRows(), m.SizeCols()) {
  // *this = m; }

  // Move constructor from Matrix
  Matrix(Matrix&& m) : MatrixView<T, ORD>(0, 0, nullptr) {
    std::swap(cols_, m.cols_);
    std::swap(rows_, m.rows_);
    std::swap(data_, m.data_);
    std::swap(d_r_, m.d_r_);
    std::swap(d_c_, m.d_c_);
  }

  template <typename TB>
  Matrix(const MatExpr<TB>& m) : Matrix(m.SizeRows(), m.SizeCols()) {
    *this = m;
  }

  ~Matrix() { delete[] data_; }

  // Copy assignment operator from Matrix
  using BASE::operator=;
  Matrix& operator=(const Matrix& m2) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        if constexpr (ORD == RowMajor) {
          data_[j + cols_ * i] = m2(i, j);
        } else {
          data_[j * rows_ + i] = m2(i, j);
        }
      }
    }

    return *this;
  }

  // Copy move assignment operator from Matrix
  Matrix& operator=(Matrix&& m2) {
    for (size_t i = 0; i < this->SizeRows(); i++) {
      for (size_t j = 0; j < this->SizeCols(); j++) {
        {
          // print(i, j)
          if constexpr (ORD == RowMajor) {
            data_[cols_ * i + j] = m2(i, j);
          } else {
            data_[i + rows_ * j] = m2(i, j);
          }
          // data_[i + cols_ * j] = m2(i, j);
        }
      }
    }

    return *this;
  };
};

template <typename... Args>
std::ostream& operator<<(std::ostream& ost, const MatrixView<Args...>& m) {
  ost << std::endl;
  for (size_t i = 0; i < m.SizeRows(); i++) {
    for (size_t j = 0; j < m.SizeCols(); j++) {
      ost << m(i, j) << " ";
    }
    ost << std::endl;
  }
  return ost;
}

constexpr ORDERING operator!(ORDERING ordering) {
  return ordering == ColMajor ? RowMajor : ColMajor;
}

template <typename T, ORDERING ORD>
auto Transpose(const Matrix<T, ORD>& m) {
  return MatrixView<T, !ORD>(m.SizeCols(), m.SizeRows(), m.Data());
}

template <typename T, ORDERING ORD>
auto Transpose(const MatrixView<T, ORD>& m) {
  return MatrixView<T, !ORD>(m.SizeCols(), m.SizeRows(), m.DistRows(),
                             m.DistCols(), m.Data());
}

template <typename T, ORDERING ORD>
auto InverseII(const Matrix<T, ORD>& m) {  // TODO: implement for the Transpose
  size_t L = m.SizeCols();
  Matrix<T, ORD> eye(L, L);
  Matrix<T, ORD> A(m);

  eye = 0;

  for (size_t i = 0; i < m.SizeRows(); i++) eye(i, i) = 1;

  for (size_t i = 0; i < eye.SizeRows(); i++) {
    T pivot = A(i, i);

    for (size_t j = 0; j < L; j++) {
      eye(i, j) = eye(i, j) / pivot;
      A(i, j) = A(i, j) / pivot;
    }

    for (size_t k = 0; k < L; k++) {
      if (k != i) {
        T factor = A(k, i);

        for (size_t j = 0; j < L; j++) {
          eye(k, j) = eye(k, j) - factor * eye(i, j);
          A(k, j) = A(k, j) - factor * A(i, j);
        }
      }
    }
  }
  return Transpose(eye);
}

template <typename T, ORDERING ORD>
auto Inverse(const Matrix<T, ORD>& m) {  // TODO: implement for the Transpose
  size_t L = m.SizeCols();
  Matrix<T, ORD> eye(L, L);
  Matrix<T, ORD> A(m);

  if constexpr (std::is_same<T, dcomplex>::value) {
    eye = dcomplex(0, 0);
    eye.Diag() = dcomplex(1, 0);
  } else {
    eye = 0;
    eye.Diag() = 1;
  }

  for (size_t i = 0; i < L; i++) {
    T pivot = A(i, i);

    // std::cout << "Row/=" << std::endl;
    eye.Row(i) /= pivot;
    A.Row(i) /= pivot;

    for (size_t j = 0; j < L; j++) {
      if (j != i) {
        T factor = A(j, i);

        eye.Row(j) = eye.Row(j) - factor * eye.Row(i);
        A.Row(j) = A.Row(j) - factor * A.Row(i);
      }
    }
  }
  return eye;
}

// in the class VectorView<T, TDIST> we have a AsMatrix(size_t rows, size_t
// cols) const method declared, we need to implement it
template <typename T, typename TDIST>
auto VectorView<T, TDIST>::AsMatrix(size_t rows, size_t cols) const {
  return MatrixView<T>(rows, cols, dist_, 1, data_);

}  // namespace Tombino_bla

// for matrix
template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
auto InnerProduct(const Matrix<TA, ORDA>& a, const Matrix<TB, ORDB>& b) {
  typedef decltype(std::declval<TA>() * std::declval<TB>()) TRES;
  TRES sum = a(0, 0) * b(0, 0);
  for (size_t i = 0; i < a.SizeRows(); i++) {
    for (size_t j = 0; j < a.SizeCols(); j++) {
      sum = sum + (a(i, j) * b(i, j));
    }
  }
  return sum;
}
// matrix view
template <typename TA, typename TB, ORDERING ORDA, ORDERING ORDB>
auto InnerProduct(const MatrixView<TA, ORDA>& a,
                  const MatrixView<TB, ORDB>& b) {
  typedef decltype(std::declval<TA>() * std::declval<TB>()) TRES;
  TRES sum = a(0, 0) * b(0, 0);
  for (size_t i = 0; i < a.SizeRows(); i++) {
    for (size_t j = 0; j < a.SizeCols(); j++) {
      sum = sum + (a(i, j) * b(i, j));
    }
  }
  return sum;
}

}  // namespace Tombino_bla
#endif
