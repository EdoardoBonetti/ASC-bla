#include <iostream>

#include <vector.h>

namespace bla = ASC_bla;


int main()
{
  size_t n = 10;
  bla::Vector<double> x(n), y(n);

  for (size_t i = 0; i < x.Size(); i++)
    {
      x(i) = i;
      y(i) = 10;
    }

  bla::Vector<double> z = x+y;
  
  std::cout << "x+y = " << z << std::endl;
  
  std::cout << "type of (x+3*y) is  " << typeid(x+3*y).name() << std::endl;

  // print the type of x+y
  x.Range(2,9) = 3;
  x.Slice(1,5) = 10;

  std::cout << "x = " << x << std::endl;

  bla::Vector<double> a(4), b(4);
  for (size_t i = 0; i < a.Size(); i++) {
    a(i) = i;
    b(i) = 100;
  }
  double sum = 0;
  sum = a * b;
  std::cout << "sum = " << sum << std::endl;
}


