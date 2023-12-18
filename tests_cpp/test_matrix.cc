#include <matrix.h>
#include <vector.h>

#include <iostream>

using namespace Tombino_bla;
using namespace std;
using std::cout;
using std::endl;
// using random numbers
#include <stdlib.h> /* srand, rand */

// create a test for the sum of two vectors

// write some test to check if the matrix is correctly initialized
int main()
{
  // define a matrix
  int n = 3;
  Matrix<double, RowMajor> A(n, n);
  Matrix<double, RowMajor> invA(n, n);

  // -4 -6 2
  // -5 -1 3
  //  -2 4 -3
  A(0, 0) = -4;
  A(0, 1) = -6;
  A(0, 2) = 2;
  A(1, 0) = -5;
  A(1, 1) = -1;
  A(1, 2) = 3;
  A(2, 0) = -2;
  A(2, 1) = 4;
  A(2, 2) = -3;

  // print it
  cout << "A = " << endl;
  cout << A << endl;

  // compute the inverse
  cout << "A^-1 = " << endl;
  invA = Inverse(A);
  cout << invA << endl;

  //  S = {{-9/118, -5/59, -8/59}, {-21/118, 8/59, 1/59}, {-11/59, 14/59,
  //  -13/59}}
  Matrix<double, RowMajor> S(n, n);
  S(0, 0) = -9.0 / 118.0;
  S(0, 1) = -5.0 / 59.0;
  S(0, 2) = -8.0 / 59.0;
  S(1, 0) = -21.0 / 118.0;
  S(1, 1) = 8.0 / 59.0;
  S(1, 2) = 1.0 / 59.0;
  S(2, 0) = -11.0 / 59.0;
  S(2, 1) = 14.0 / 59.0;
  S(2, 2) = -13.0 / 59.0;

  // print S- A^-1
  cout << "S  = " << endl;
  cout << S << endl;

  cout << "S - A^-1 = " << endl;
  cout << S - invA << endl;

  // check if the entries are almost zero
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; ++j) {
      if (abs(S(i, j) - invA(i, j)) > 1e-10) {
        cout << "Error: the entries are not almost zero" << endl;
        return 1;
      }
    }
  }

  return 0;
}