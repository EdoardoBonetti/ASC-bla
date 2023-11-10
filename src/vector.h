#ifndef FILE_VECTOR_H
#define FILE_VECTOR_H

#include <iostream>

#include "expression.h"

namespace Tombino_bla {

template <typename T, typename TDIST = std::integral_constant<size_t, 1>>
class VectorView : public VecExpr<VectorView<T, TDIST>> {
 protected:
  T* data_;
  size_t size_;
  TDIST dist_;

 public:
  VectorView(size_t size, T* data) : data_(data), size_(size) {}

  VectorView(size_t size, TDIST dist, T* data) : data_(data), size_(size), dist_(dist) {}

  template <typename TB>
  VectorView& operator=(const VecExpr<TB>& v2) {
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

  auto Range(size_t first, size_t next) const { return VectorView(next - first, dist_, data_ + first); }

  auto Slice(size_t first, size_t slice) const {
    return VectorView<T, size_t>(size_ / slice, dist_ * slice, data_ + first * dist_);
  }
};

  template <typename T>
  class Vector : public VectorView<T>
  {
    typedef VectorView<T> BASE;
    using BASE::size_;
    using BASE::data_;
  public:
    Vector (size_t size) 
      : VectorView<T> (size, new T[size]) { ; }
    
    Vector (const Vector & v)
      : Vector(v.Size())
    {
      *this = v;
    }

    Vector (Vector && v)
      : VectorView<T> (0, nullptr)
    {
      std::swap(size_, v.size_);
      std::swap(data_, v.data_);
    }

    template <typename TB>
    Vector (const VecExpr<TB> & v)
      : Vector(v.Size())
    {
      *this = v;
    }
    
    
    ~Vector () { delete [] data_; }

    using BASE::operator=;
    Vector & operator=(const Vector & v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }

    Vector & operator= (Vector && v2)
    {
      for (size_t i = 0; i < size_; i++)
        data_[i] = v2(i);
      return *this;
    }
    
  };


  template <typename ...Args>
  std::ostream & operator<< (std::ostream & ost, const VectorView<Args...> & v)
  {
    if (v.Size() > 0)
      ost << v(0);
    for (size_t i = 1; i < v.Size(); i++)
      ost << ", " << v(i);
    return ost;
  }

  }  // namespace Tombino_bla

#endif
