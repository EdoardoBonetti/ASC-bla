#include <vector.h>

#include <cassert>
#include <iostream>
#include <vector>

using namespace Tombino_bla;
using namespace std;

// create a test for the sum of two vecs
auto test_vec_sum()
{
  Vec<3> v1(3);
  cout << v1 << endl;
  Vec<3> v2(3);

  return 0;
}

int main()
{
  // run the test
  assert(test_vec_sum() == 0);
  return 0;
}