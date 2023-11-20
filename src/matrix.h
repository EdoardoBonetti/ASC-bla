#ifndef FILE_MATRIX_H
#define FILE_MATRIX_H

#include <iostream>

#include "expression.h"
#include "vector.h"

namespace Tombino_bla {
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
    if constexpr (ORD == ColMajor) {
      dist_ = cols;  // rows;
    } else {
      dist_ = rows;
    }
  }

  MatrixView(size_t rows, size_t cols, size_t dist, T* data)
      : data_(data), rows_(rows), cols_(cols), dist_(dist) {}

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
        std::cout << "i = " << i << " j = " << j << std::endl;
        (*this)(i, j) = scal;
      }
    }
    return *this;
  }

  auto View() const { return MatrixView(rows_, cols_, dist_, data_); }
  size_t SizeCols() const { return cols_; }
  size_t SizeRows() const { return rows_; }
  T* Data() const { return data_; }
  size_t Dist() const { return dist_; }
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
    if constexpr (ORD == RowMajor) {
      return VectorView<T>(rows_, data_ + i * dist_);
    } else {
      return VectorView<T, size_t>(rows_, dist_, data_ + i);
    }
  }

  auto Rows(size_t first, size_t next) const {
    if constexpr (ORD == RowMajor) {
      return MatrixView<T, ORD>(next - first, cols_, dist_, data_ + first * dist_);
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

  // Copy constructor from Matrix
  Matrix(const Matrix& m) : Matrix(m.SizeRows(), m.SizeCols()) { *this = m; }

  // Copy constructor from MatrixView: New feature
  // Matrix(const MatrixView<T, ORD>& m) : Matrix(m.SizeRows(), m.SizeCols()) { *this = m; }

  // Move constructor from Matrix
  Matrix(Matrix&& m) : MatrixView<T, ORD>(0, 0, nullptr) {
    std::swap(cols_, m.cols_);
    std::swap(rows_, m.rows_);
    std::swap(data_, m.data_);
    std::swap(dist_, m.dist_);
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
        std::cout << "i = " << i << " j = " << j << std::endl;

        if constexpr (ORD == RowMajor) {
          data_[i + cols_ * j] = m2(i, j);
        } else {
          data_[i * rows_ + j] = m2(i, j);
        }
        // data_[i + cols_ * j] = m2(i, j);
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
          std::cout << "i = " << i << " j = " << j << std::endl;
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
auto Transpose(const MatrixView<T, ORD>& m) {
  return MatrixView<T, !ORD>(m.SizeCols(), m.SizeRows(), m.Dist(), m.Data());
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

}  // namespace Tombino_bla
#endif
