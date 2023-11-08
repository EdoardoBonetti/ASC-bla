#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"
#include "vector.h"

namespace ASC_bla {
enum ORDERING { RowMajor, ColMajor };
template <typename T, ORDERING ORD>
class MatrixView : public MatExpr<MatrixView<T, ORD>> {
 protected:
  T* data_;
  size_t cols_;
  size_t rows_;
  size_t dist_;

 public:
  // Constructor
  MatrixView(size_t rows, size_t cols, T* data)
      : data_(data), rows_(rows), cols_(cols) {
    if (ORD == ColMajor) {
      dist_ = cols;
    } else {
      dist_ = rows;
    }
  }

  MatrixView(size_t rows, size_t cols, size_t dist, T* data)
      : data_(data), rows_(rows), cols_(cols), dist_(dist) {}

  MatrixView& operator=(const MatrixView& m2) {
    *this = static_cast<const MatExpr<MatrixView<T, ORD>>&>(m2);
    return *this;
  }
  // Copy assignment operator, needs to be modified for row maj
  template <typename TB>
  MatrixView& operator=(const MatExpr<TB>& m2) {
    for (size_t i = 0; i < rows_; i++) {
      for (size_t j = 0; j < cols_; j++) {
        (*this)(i, j) = m2(i, j);
      }
    }
    return *this;
  }

  MatrixView& operator=(T scal) {
    for (size_t i = 0; i < rows_; i++) {
      for (size_t j = 0; j < cols_; j++) {
        (*this)(i, j) = scal;
      }
    }
    return *this;
  }

  auto View() const { return MatrixView(rows_, cols_, dist_, data_); }
  size_t SizeCols() const { return cols_; }
  size_t SizeRows() const { return rows_; }
  T* Data() const { return data_; }
  T& operator()(size_t i, size_t j) {
    if (ORD == RowMajor) {
      return data_[dist_ * i + j];
    } else {
      return data_[dist_ * j + i];
    }
  }
  const T& operator()(size_t i, size_t j) const {
    if (ORD == RowMajor) {
      return data_[dist_ * i + j];
    } else {
      return data_[dist_ * j + i];
    }
  }

  auto Row(size_t i) const {
    if constexpr (ORD == RowMajor) {
      return VectorView<T>(cols_, data_ + i * dist_);
    } else {
      return VectorView<T, size_t>(cols_, dist_, data_ + i);
    }
  }

  auto Col(size_t i) const {
    if constexpr (ORD == ColMajor) {
      return VectorView<T>(rows_, data_ + i * dist_);
    } else {
      return VectorView<T, size_t>(rows_, dist_, data_ + i);
    }
  }

  auto Rows(size_t first, size_t next) const {
    if constexpr (ORD == ColMajor) {
      return MatrixView<T, ORD>(next - first, cols_, dist_,
                                data_ + first * dist_);
    } else {
      return MatrixView<T, ORD>(next - first, cols_, dist_, data_ + first);
    }
  }

  auto Cols(size_t first, size_t next) const {
    if constexpr (ORD == RowMajor) {
      return MatrixView<T, ORD>(rows_, next - first, dist_, data_ + first);
    } else {
      return MatrixView<T, ORD>(rows_, next - first, dist_,
                                data_ + first * dist_);
    }
  }
};

template <typename T, ORDERING ORD>
class Matrix : public MatrixView<T, ORD> {
  typedef MatrixView<T, ORD> BASE;
  using BASE::cols_;
  using BASE::data_;
  using BASE::dist_;
  using BASE::rows_;

 public:
  Matrix(size_t rows, size_t cols)
      : MatrixView<T, ORD>(rows, cols, new T[rows * cols]) {
    ;
  }

  Matrix(const Matrix& m) : Matrix(m.SizeRows(), m.SizeCols()) { *this = m; }

  Matrix(Matrix&& m) : MatrixView<T, ORD>(0, 0, nullptr) {
    std::swap(cols_, m.cols_);
    std::swap(rows_, m.rows_);
    std::swap(data_, m.data_);
    std::swap(dist_, m.dist_);
  }

  template <typename TB>
  Matrix(const MatExpr<TB>& m) : Matrix(m.SizeCols(), m.SizeRows()) {
    *this = m;
  }

  ~Matrix() { delete[] data_; }

  using BASE::operator=;
  Matrix& operator=(const Matrix& m2) {
    for (size_t i = 0; i < cols_; i++) {
      for (size_t j = 0; j < rows_; j++) {
        {
          data_[i + cols_ * j] = m2(i, j);
        }
      }
    }

    return *this;
  }

  Matrix& operator=(Matrix&& m2) {
    for (size_t i = 0; i < cols_; i++) {
      for (size_t j = 0; j < rows_; j++) {
        {
          data_[i + cols_ * j] = m2(i, j);
        }
      }
    }

    return *this;
  };
};

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

constexpr ORDERING operator!(ORDERING ordering) {
  return ordering == ColMajor ? RowMajor : ColMajor;
}

template <typename T, ORDERING ORD>
auto Transpose(const Matrix<T, ORD>& m) {
  return MatrixView<T, !ORD>(m.SizeCols(), m.SizeRows(), m.Data());
}

template <typename T, ORDERING ORD>
auto Inverse(const Matrix<T, ORD>& m) {
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
  return eye;
}

}  // namespace ASC_bla
#endif
