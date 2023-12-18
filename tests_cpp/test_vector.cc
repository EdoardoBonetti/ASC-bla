#include <vector.h>

#include <cassert>
#include <iostream>
#include <vector>

using namespace Tombino_bla;
using namespace std;

// create a test for the sum of two vectors
auto test_vector_sum()
{
  // Your Vector class operations and assertions here
  Vector<double> v1(3);
  Vector<double> v2(3);
  v1(0) = 1;
  v1(1) = 2;
  v1(2) = 3;
  v2(0) = 4;
  v2(1) = 5;
  v2(2) = 6;
  Vector result = v1 + v2;
  if (result(0) == 5 && result(1) == 7 && result(2) == 9)
    return 0;
  else
    return 1;
}

// write some test to check if the matrix is correctly initialized
int main()
{
  // run the test
  assert(test_vector_sum() == 0);
  return 0;
}