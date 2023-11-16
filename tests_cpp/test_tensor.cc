#include <iostream>

#include "tensor.h"

namespace bla = Tombino_bla;

int main() {
  // Example usage
  Tensor<int, 3, 4, 5> myTensor;

  // Set values using the operator (i, j) syntax
  myTensor(0, 0, 0) = 1;
  myTensor(1, 2, 2) = 42;

  // Print the tensor
  myTensor.print();

  return 0;
}
