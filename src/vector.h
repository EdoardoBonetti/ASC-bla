#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>

#include "expression.h"
// #include "matrix.h"

namespace Tombino_bla {

auto Complement(dcomplex& v) { return v; }

template <typename T = double,
          typename TDIST = std::integral_constant<size_t, 1>>
class VectorView : public VecExpr<VectorView<T, TDIST>> {
 protected:
  T* data_;
  size_t size_;
  TDIST dist_;

 public:
  VectorView(size_t size) : data_(new T[size]), size_(size) {}  // new
  VectorView(size_t size, T* data) : data_(data), size_(size) {}
  VectorView(size_t size, TDIST dist)
      : data_(new T[size * dist]), size_(size), dist_(dist) {}  // new
  VectorView(size_t size, TDIST dist, T* data)
      : data_(data), size_(size), dist_(dist) {}
  VectorView(std::initializer_list<T> list)
      : data_(new T[list.size()]), size_(list.size()) {
    std::copy(list.begin(), list.end(), data_);
  }

  template <typename TB>
  VectorView& operator=(const VecExpr<TB>& v2) {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = v2(i);
    return *this;
  }

  VectorView& operator=(const VectorView& v2) {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = v2(i);
    return *this;
  }

  VectorView& operator=(T scal) {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = scal;
    return *this;
  }

  auto View() const { return VectorView(size_, dist_, data_); }
  size_t Size() const { return size_; }
  size_t Dist() const { return dist_; }
  T& operator()(size_t i) { return data_[dist_ * i]; }
  const T& operator()(size_t i) const { return data_[dist_ * i]; }

  // Data()
  T* Data() { return data_; }
  const T* Data() const { return data_; }

  auto Range(size_t first, size_t next) const {
    return VectorView(next - first, dist_, data_ + first * dist_);
  }
  auto Slice(size_t first, size_t slice) const {
    return VectorView<T, size_t>((size_ + 1 - first) / slice, dist_ * slice,
                                 data_ + first * dist_);
  }

  // operator+= vec
  template <typename TB>
  VectorView& operator+=(const VecExpr<TB>& v2) {
    for (size_t i = 0; i < size_; i++) this->operator()(i) += v2(i);
    return *this;
  }

  // operator-= vec
  template <typename TB>
  VectorView& operator-=(const VecExpr<TB>& v2) {
    for (size_t i = 0; i < size_; i++) this->operator()(i) -= v2(i);
    return *this;
  }

  // operator*= scal
  template <ValidSCAL TSCAL>
  VectorView& operator*=(const TSCAL& scal) {
    for (size_t i = 0; i < size_; i++) this->operator()(i) *= scal;
    return *this;
  }

  // operator/= scal
  template <ValidSCAL TSCAL>
  VectorView& operator/=(const TSCAL& scal) {
    for (size_t i = 0; i < size_; i++) this->operator()(i) /= scal;
    return *this;
  }

  // operator==
  template <typename TB, typename TDIST2>
  bool operator==(const VectorView<TB, TDIST2>& v2) const {
    for (size_t i = 0; i < size_; i++)
      if (this->operator()(i) != v2(i)) {
        // get rid of next 3 lines
        std::cout << "operator== failed" << std::endl;
        std::cout << v2(i) << std::endl;
        std::cout << this->operator()(i) << std::endl;
        return false;
      }
    return true;
  }

  // operator!=
  template <typename TB, typename TDIST2>
  bool operator!=(const VectorView<TB, TDIST2>& v2) const {
    for (size_t i = 0; i < size_; i++)
      if (this->operator()(i) != v2(i)) {
        // get rid of next 3 lines
        std::cout << "operator!= failed" << std::endl;
        std::cout << v2(i) << std::endl;
        std::cout << this->operator()(i) << std::endl;
        return true;
      }

    return false;
  }

  // AsMatrix
  auto AsMatrix(size_t rows, size_t cols) const;
};

template <typename T = double>
class Vector : public VectorView<T> {
  typedef VectorView<T> BASE;
  using BASE::data_;
  using BASE::size_;

 public:
  Vector(size_t size) : VectorView<T>(size, new T[size]) { ; }
  Vector(size_t size, T* data) : VectorView<T>(size, data) { ; }

  Vector(const Vector& v) : Vector(v.Size()) { *this = v; }

  // list initialization
  Vector(std::initializer_list<T> list) : VectorView<T>(list) {}

  Vector(Vector&& v) : VectorView<T>(0, nullptr) {
    std::swap(size_, v.size_);
    std::swap(data_, v.data_);
  }

  template <typename TB>
  Vector(const VecExpr<TB>& v) : Vector(v.Size()) {
    *this = v;
  }

  ~Vector() { delete[] data_; }

  using BASE::operator=;
  Vector& operator=(const Vector& v2) {
    for (size_t i = 0; i < size_; i++) data_[i] = v2(i);
    return *this;
  }

  Vector& operator=(Vector&& v2) {
    for (size_t i = 0; i < size_; i++) data_[i] = v2(i);
    return *this;
  }
};

template <int SIZE = 3, typename T = double>
class Vec : public VecExpr<Vec<SIZE, T>> {
 private:
  T data_[SIZE];

 public:
  Vec() = default;

  Vec(const Vec& v2) {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }

  template <typename TB>
  Vec(const VectorView<TB>& v2) {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }
  template <typename TB>
  Vec(const VecExpr<TB>& v2) {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }

  Vec(std::initializer_list<T> list) {
    for (size_t i = 0; i < list.size(); i++) data_[i] = list.begin()[i];
  }

  Vec(T scal) {
    for (size_t i = 0; i < SIZE; i++) data_[i] = scal;
  }

  auto Upcast() const { return *this; }
  size_t Size() const { return SIZE; }
  T& operator()(size_t i) { return data_[i]; }
  const T& operator()(size_t i) const { return data_[i]; }
};

template <typename... Args>
std::ostream& operator<<(std::ostream& ost, const VectorView<Args...>& v) {
  if (v.Size() > 0) ost << v(0);
  for (size_t i = 1; i < v.Size(); i++) ost << ", " << v(i);
  return ost;
}

// define the InnerProduct for Vector and VectorView
template <typename T>
auto InnerProduct(const VectorView<T>& v1, const VectorView<T>& v2) {
  T sum = 0;
  for (size_t i = 0; i < v1.Size(); i++) {
    if constexpr (std::is_same<T, dcomplex>::value)
      sum += v1(i) * Complement(v2(i));
    else
      sum += v1(i) * v2(i);
  }

  return sum;
}

// inner product for VectorView
template <typename TA, typename TB>
auto InnerProduct(const VectorView<TA>& a, const VectorView<TB>& b) {
  typedef decltype(std::declval<TA>() * std::declval<TB>()) TRES;
  TRES sum = a(0) * b(0);
  for (size_t i = 1; i < a.Size(); i++) sum = sum + (a(i) * b(i));
  return sum;
}
// for Vector
template <typename TA, typename TB>
auto InnerProduct(const Vector<TA>& a, const Vector<TB>& b) {
  typedef decltype(std::declval<TA>() * std::declval<TB>()) TRES;
  TRES sum = a(0) * b(0);
  for (size_t i = 1; i < a.Size(); i++) sum = sum + (a(i) * b(i));
  return sum;
}

}  // namespace Tombino_bla

#endif
