#include <iostream>
#include <vector>

namespace Tombino_bla {

template <typename... Indices>
size_t index(Indices... indices) const {
  static_assert(sizeof...(Indices) == sizeof...(Dims), "Number of indices must match number of dimensions");
  size_t indicesArray[] = {static_cast<size_t>(indices)...};
  size_t linearIndex = 0;
  size_t multiplier = 1;
  size_t dimSizes[] = {Dims...};
  for (size_t i = 0; i < sizeof...(Dims); ++i) {
    linearIndex += indicesArray[i] * multiplier;
    multiplier *= dimSizes[i];
  }
  return linearIndex;
}

template <typename T, typename TDIST = std::integral_constant<size_t, 0>>
class TensorShape {
  // tensor shape is an auxiliary class that stores the shape of a tensor:
  // - number of dimensions
  // - size of each dimension

 private:
  size_t ndims_;
  size_t* dims_;

 public:
  TensorShape(size_t ndims) : ndims_(ndims) { dims_ = new size_t[ndims_]; }
  TensorShape(size_t* dims) : ndims_(sizeof(dims) / sizeof(dims[0])), dims_(dims) {}
  TensorShape(size_t ndims, size_t* dims) : ndims_(ndims), dims_(dims) {}

  size_t ndims() const { return ndims_; }
  size_t* dims() const { return dims_; }
};

template <typename T>
class TensorView {
 private:
  T* data;
  TesnsorShape shape;
  TensorShape dist;

 public:
  TensorView(size_t ndims, size_t* dims, T* data) : shape(ndims, dims), data_(data) {
    // initialize dist with only 1s
    size_t* dist = new size_t[ndims];
    for (size_t i = 0; i < ndims; i++) {
      dist[i] = 1;
    }
  }
  TensorView(size_t ndims, size_t* dims, TensorShape dist, T* data) : shape(ndims, dims), dist_(dist), data_(data) {}

  TensorView& operator=(const TensorView& t2) {
    *this = static_cast<const TensExpr<TensorView<T, ORD>>&>(t2);
    return *this;
  }

  // TensorView& operator=(const T scal) {
  //   for (size_t i = 0; i < size_; i++) {
  //     data_[dist_ * i] = scal;
  //   }
  //   return *this;
  // }

  // first we define the operaror() using using the operator (i, j, ...) syntax but knowing tht there are ndims indices

  // Tensor(
  // Calculate the total size of the tensor based on the dimensions
  //     size_t size = (Dims * ...);
  //     data.resize(size);
  //   }
  //
  //   // Helper function to convert multi-dimensional indices to a linear index
  //   template <typename... Indices>
  //   size_t index(Indices... indices) const {
  //     static_assert(sizeof...(Indices) == sizeof...(Dims), "Number of indices must match number of dimensions");
  //     size_t indicesArray[] = {static_cast<size_t>(indices)...};
  //     size_t linearIndex = 0;
  //     size_t multiplier = 1;
  //     size_t dimSizes[] = {Dims...};
  //     for (size_t i = 0; i < sizeof...(Dims); ++i) {
  //       linearIndex += indicesArray[i] * multiplier;
  //       multiplier *= dimSizes[i];
  //     }
  //     return linearIndex;
  //   }
  //
  //   // Accessor using the operator (i, j, ...) syntax
  //   template <typename... Indices>
  //   T &operator()(Indices... indices) {
  //     return data[index(indices...)];
  //   }
  //
  //   // Const accessor
  //   template <typename... Indices>
  //   const T &operator()(Indices... indices) const {
  //     return data[index(indices...)];
  //   }
  //
  //   // Example function to print the tensor
  //   void print() const {
  //     for (const auto &value : data) {
  //       std::cout << value << " ";
  //     }
  //     std::cout << std::endl;
  //   }
};

}  // namespace Tombino_bla

int main() {
  //// Example usage
  // Tensor<int, 3, 4, 5> myTensor;
  //
  //// Set values using the operator (i, j) syntax
  // myTensor(0, 0, 0) = 1;
  // myTensor(1, 2, 2) = 42;
  //
  //// Print the tensor
  // myTensor.print();
  //
  TensorShape shape();

  return 0;
}
