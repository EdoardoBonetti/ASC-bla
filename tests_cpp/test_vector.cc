// File to test the Vector.h file
//
// Author: Edoardo Bonetti
//
// Date: 20/12/2023

#include <vector.h>

#include <cassert>
#include <iostream>
#include <vector>

using namespace Tombino_bla;
using namespace std;

// the outline of a test must be the following:
/*
auto test_<name_object/s>_<name_operation>()
{

  create the object/s
  perform the operation
  assert the result against the expected one

    - return 1 if fails assert and print result
    - repeat thea above for all possibilities

  return 0 if all tests pass

}
*/

/********************************************/
/*     TESTS FOR THE VECTOR VIEW CLASS      */
/********************************************/

// test vector view constructor
auto test_vectorview_constructor() {
  // Constructor
  VectorView<double> v(3);
  v(0) = 1;
  v(1) = 2;
  v(2) = 3;

  if (v.Size() != 3 || v.Dist() != 1) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Size = " << v.Size() << endl;
    cout << "Dist = " << v.Dist() << endl;
    return 1;
  }

  if (v.Data() != &v(0) || v.Data() + 1 != &v(1) || v.Data() + 2 != &v(2)) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Data = " << v.Data() << endl;
    cout << "Data()+1 = " << v.Data() + 1 << endl;
    cout << "Data()+2 = " << v.Data() + 2 << endl;
    return 1;
  }

  // constructor with data
  double* data = new double[3];
  data[0] = 1;
  data[1] = 2;
  data[2] = 3;
  VectorView<double> v1(3, data);

  if (v1.Size() != 3 || v1.Dist() != 1) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Size = " << v.Size() << endl;
    cout << "Dist = " << v.Dist() << endl;
    return 1;
  }

  if (v1.Data() != &v1(0) || v1.Data() + 1 != &v1(1) ||
      v1.Data() + 2 != &v1(2)) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Data = " << v1.Data() << endl;
    cout << "Data()+1 = " << v1.Data() + 1 << endl;
    cout << "Data()+2 = " << v1.Data() + 2 << endl;
    return 1;
  }

  // constructor with data and dist = 2
  double* data2 = new double[3];
  data2[0] = 1;
  data2[1] = 2;
  data2[2] = 3;
  size_t dist = 2;
  size_t size = 3;
  VectorView<double, size_t> v2(size, dist, data2);

  if (v2.Size() != 3 || v2.Dist() != 2) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Size = " << v2.Size() << endl;
    cout << "Dist = " << v2.Dist() << endl;
    return 1;
  }

  if (v2.Data() != &v2(0) || v2.Data() + 2 != &v2(1) ||
      v2.Data() + 4 != &v2(2)) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Data = " << v2.Data() << endl;
    cout << "Data()+1 = " << v2.Data() + 2 << endl;
    cout << "Data()+2 = " << v2.Data() + 4 << endl;
    return 1;
  }

  // constructor with data and dist = 2
  double* data3 = new double[3];
  data3[0] = 1;
  data3[1] = 2;
  data3[2] = 3;
  size_t dist3 = 3;
  size_t size3 = 3;
  VectorView<double, size_t> v3(size3, dist3, data3);

  if (v3.Size() != 3 || v3.Dist() != 3) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Size = " << v3.Size() << endl;
    cout << "Dist = " << v3.Dist() << endl;
    return 1;
  }

  if (v3.Data() != &v3(0) || v3.Data() + 3 != &v3(1) ||
      v3.Data() + 6 != &v3(2)) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Data = " << v3.Data() << endl;
    cout << "Data()+1 = " << v3.Data() + 3 << endl;
    cout << "Data()+2 = " << v3.Data() + 6 << endl;
    return 1;
  }

  // cons

  return 0;
}

auto test_vectorview_sum() {
  // Sum of vectors
  VectorView<double> v1(3);
  VectorView<double> v2(3);
  double alpha = 2.0;
  double beta = 3.0;
  v1(0) = 1;
  v1(1) = 2;
  v1(2) = 3;
  v2(0) = 4;
  v2(1) = 5;
  v2(2) = 6;
  Vector result = alpha * v1 + beta * v2;
  if (result(0) != 14 && result(1) != 19 && result(2) != 24) return 1;

  size_t dist = 5;
  size_t size = 3;

  VectorView<double, size_t> v3(size, dist);
  // cout << v3 << endl;

  // v3(0) = 1;
  // v3(1) = 2;
  // v3(2) = 3;

  // for (size_t i = 0; i < 15; i++) {
  //   cout << v3.Data()[i] << endl;
  // }
  //
  if (v3.Size() != 3 || v3.Dist() != 5) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Size = " << v3.Size() << endl;
    cout << "Dist = " << v3.Dist() << endl;
    return 1;
  }

  if (v3.Data() != &v3(0) || v3.Data() + 5 != &v3(1) ||
      v3.Data() + 10 != &v3(2)) {
    cout << "Error in VectorView constructor" << endl;
    cout << "Data = " << v3.Data() << endl;
    cout << "Data()+1 = " << v3.Data() + 5 << endl;
    cout << "Data()+2 = " << v3.Data() + 10 << endl;
    return 1;
  }

  return 0;
}

auto test_vectorview_range() {
  // create a vector view for every constructor
  size_t dist = 20;
  size_t size = 10;
  double* data = new double[200];

  VectorView<double> v1(size);
  VectorView<double> v2(size, data);
  VectorView<double, size_t> v3(size, dist);
  VectorView<double, size_t> v4(size, dist, data);

  // now set the entirs to 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
  for (size_t i = 0; i < size; i++) {
    v1(i) = i;
    v2(i) = i;
    v3(i) = i;
    v4(i) = i;
  }

  // test slice(0,5)
  if (v1.Range(0, 5) != v2.Range(0, 5) || v1.Range(0, 5) != v3.Range(0, 5) ||
      v1.Range(0, 5) != v4.Range(0, 5)) {
    cout << "Error in VectorView Range" << endl;
    cout << "v1.Range(0, 5) = " << v1.Range(0, 5) << endl;
    cout << "v2.Range(0, 5) = " << v2.Range(0, 5) << endl;
    cout << "v3.Range(0, 5) = " << v3.Range(0, 5) << endl;
    cout << "v4.Range(0, 5) = " << v4.Range(0, 5) << endl;
    return 1;
  }
  // test slice(1, 4)
  if (v1.Range(1, 4) != v2.Range(1, 4) || v1.Range(1, 4) != v3.Range(1, 4) ||
      v1.Range(1, 4) != v4.Range(1, 4)) {
    cout << "Error in VectorView Range" << endl;
    cout << "v1.Range(1,4) = " << v1.Range(1, 4) << endl;
    cout << "v2.Range(1,4) = " << v2.Range(1, 4) << endl;
    cout << "v3.Range(1,4) = " << v3.Range(1, 4) << endl;
    cout << "v4.Range(1,4) = " << v4.Range(1, 4) << endl;
    return 1;
  }

  // test slice(3, 7)
  if (v1.Range(3, 7) != v2.Range(3, 7) || v1.Range(3, 7) != v3.Range(3, 7) ||
      v1.Range(3, 7) != v4.Range(3, 7)) {
    cout << "Error in VectorView Range" << endl;
    cout << "v1.Range(3,7) = " << v1.Range(3, 7) << endl;
    cout << "v2.Range(3,7) = " << v2.Range(3, 7) << endl;
    cout << "v3.Range(3,7) = " << v3.Range(3, 7) << endl;
    cout << "v4.Range(3,7) = " << v4.Range(3, 7) << endl;
    return 1;
  }
  return 0;
}
/********************************************/
/*          TESTS FOR THE VECTOR CLASS      */
/********************************************/

auto test_vector_sum() {
  // Sum of vectors
  Vector<double> v1(3);
  Vector<double> v2(3);
  double alpha = 2.0;
  double beta = 3.0;
  v1(0) = 1;
  v1(1) = 2;
  v1(2) = 3;
  v2(0) = 4;
  v2(1) = 5;
  v2(2) = 6;
  Vector result = alpha * v1 + beta * v2;
  if (result(0) == 14 && result(1) == 19 && result(2) == 24)
    return 0;
  else
    return 1;
}

auto test_vector_range() { return 0; }

// write some test to check if the matrix is correctly initialized
int main() {
  // run the test
  assert(test_vectorview_constructor() == 0);
  assert(test_vectorview_sum() == 0);
  assert(test_vectorview_range() == 0);
  return 0;
}