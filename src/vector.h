#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>

#include "expression.h"
// #include "matrix.h"

namespace Tombino_bla
{

template <typename T = double,
          typename TDIST = std::integral_constant<size_t, 1>>
class VectorView : public VecExpr<VectorView<T, TDIST>>
{
 protected:
  T* data_;
  size_t size_;
  TDIST dist_;

 public:
  VectorView(size_t size, T* data) : data_(data), size_(size) {}
  VectorView(size_t size, TDIST dist, T* data)
      : data_(data), size_(size), dist_(dist)
  {
  }
  // copy constructor for const size and dist etc

  // VectorView(const size_t size, const T* data) : data_(data), size_(size) {}
  // VectorView(const size_t size, const TDIST dist, const T* data)
  //     : data_(data), size_(size), dist_(dist)
  //{
  // }

  // brace-enclosed initializer list
  VectorView(std::initializer_list<T> list)
      : data_(new T[list.size()]), size_(list.size())
  {
    std::copy(list.begin(), list.end(), data_);
  }

  template <typename TB>
  VectorView& operator=(const VecExpr<TB>& v2)
  {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = v2(i);
    return *this;
  }

  VectorView& operator=(const VectorView& v2)
  {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = v2(i);
    return *this;
  }

  VectorView& operator=(T scal)
  {
    for (size_t i = 0; i < size_; i++) data_[dist_ * i] = scal;
    return *this;
  }

  auto View() const { return VectorView(size_, dist_, data_); }
  size_t Size() const { return size_; }
  size_t Dist() const { return dist_; }
  T& operator()(size_t i) { return data_[dist_ * i]; }
  const T& operator()(size_t i) const { return data_[dist_ * i]; }

  auto Range(size_t first, size_t next) const
  {
    return VectorView(next - first, dist_, data_ + first * dist_);
  }
  auto Slice(size_t first, size_t slice) const
  {
    return VectorView<T, size_t>(size_ / slice, dist_ * slice,
                                 data_ + first * dist_);
  }

  // operator+= vec
  template <typename TB>
  VectorView& operator+=(const VecExpr<TB>& v2)
  {
    for (size_t i = 0; i < size_; i++) this->operator()(i) += v2(i);
    return *this;
  }

  // operator-= vec
  template <typename TB>
  VectorView& operator-=(const VecExpr<TB>& v2)
  {
    for (size_t i = 0; i < size_; i++) this->operator()(i) -= v2(i);
    return *this;
  }

  // operator*= scal
  template <ValidSCAL TSCAL>
  VectorView& operator*=(const TSCAL& scal)
  {
    for (size_t i = 0; i < size_; i++) this->operator()(i) *= scal;
    return *this;
  }

  // operator/= scal
  template <ValidSCAL TSCAL>
  VectorView& operator/=(const TSCAL& scal)
  {
    for (size_t i = 0; i < size_; i++) this->operator()(i) /= scal;
    return *this;
  }

  // AsMatrix
  auto AsMatrix(size_t rows, size_t cols) const;
};

template <typename T = double>
class Vector : public VectorView<T>
{
  typedef VectorView<T> BASE;
  using BASE::data_;
  using BASE::size_;

 public:
  Vector(size_t size) : VectorView<T>(size, new T[size]) { ; }

  Vector(const Vector& v) : Vector(v.Size()) { *this = v; }

  // list initialization
  Vector(std::initializer_list<T> list) : VectorView<T>(list) {}

  Vector(Vector&& v) : VectorView<T>(0, nullptr)
  {
    std::swap(size_, v.size_);
    std::swap(data_, v.data_);
  }

  template <typename TB>
  Vector(const VecExpr<TB>& v) : Vector(v.Size())
  {
    *this = v;
  }

  ~Vector() { delete[] data_; }

  using BASE::operator=;
  Vector& operator=(const Vector& v2)
  {
    for (size_t i = 0; i < size_; i++) data_[i] = v2(i);
    return *this;
  }

  Vector& operator=(Vector&& v2)
  {
    for (size_t i = 0; i < size_; i++) data_[i] = v2(i);
    return *this;
  }
};

template <int SIZE = 3, typename T = double>
class Vec : public VecExpr<Vec<SIZE, T>>
{
 private:
  T data_[SIZE];

 public:
  Vec() = default;
  // Vec(T scal) { *this = scal; }
  // Vec(std::initializer_list<T> list) { *this = list; }

  // copy constructor
  // Vec(const Vec& v) { *this = v; }
  // copy constructor from VectorView
  // template <typename... Args>
  // Vec(const VectorView<Args...>& v) : Vec()
  //{
  //  for (size_t i = 0; i < SIZE; i++) data[i] = v(i);
  //}
  // copy constructor from VecExpr
  // template <typename TB>
  // Vec(const VecExpr<TB>& v) : Vec()
  //{
  //  for (size_t i = 0; i < SIZE; i++) data[i] = v(i);
  //}

  // move constructor
  // Vec(Vec&& v) { *this = std::move(v); }

  // assignment operator
  // Vec& operator=(const Vec& v)
  //{
  //  for (int i = 0; i < SIZE; i++) data[i] = v(i);
  //  return *this;
  //}

  // move assignment operator
  // Vec& operator=(Vec&& v)
  //{
  //  for (int i = 0; i < SIZE; i++) data[i] = v(i);
  //  return *this;
  //}

  // access operator
  // T& operator()(size_t i) { return data[i]; }
  // const T& operator()(size_t i) const { return data[i]; }
  //
  //// operator+= vec
  // template <typename TB>
  // Vec& operator+=(const VecExpr<TB>& v2)
  //{
  //   for (size_t i = 0; i < SIZE; i++) this->operator()(i) += v2(i);
  //   return *this;
  // }

  Vec(const Vec& v2)
  {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }

  template <typename TB>
  Vec(const VectorView<TB>& v2)
  {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }
  template <typename TB>
  Vec(const VecExpr<TB>& v2)
  {
    for (size_t i = 0; i < SIZE; i++) data_[i] = v2(i);
  }

  Vec(std::initializer_list<T> list)
  {
    for (size_t i = 0; i < list.size(); i++) data_[i] = list.begin()[i];
  }

  Vec(T scal)
  {
    for (size_t i = 0; i < SIZE; i++) data_[i] = scal;
  }

  auto Upcast() const { return *this; }
  size_t Size() const { return SIZE; }
  T& operator()(size_t i) { return data_[i]; }
  const T& operator()(size_t i) const { return data_[i]; }
};

template <typename... Args>
std::ostream& operator<<(std::ostream& ost, const VectorView<Args...>& v)
{
  if (v.Size() > 0) ost << v(0);
  for (size_t i = 1; i < v.Size(); i++) ost << ", " << v(i);
  return ost;
}

}  // namespace Tombino_bla

#endif
