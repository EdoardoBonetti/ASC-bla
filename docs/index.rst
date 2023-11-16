.. ASC-bla documentation master file, created by
   sphinx-quickstart on Tue Aug 29 06:39:02 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to TomBino's documentation!
===================================

ASC-bla is a C++ library for basic linear algebra operations.
The library provides template classes **Vector** and **Matrix**.

Installation is via git-clone:

..  code-block::
    
    git clone https://github.com/EdoardoBonetti/ASC-bla


To configure and build some tests do

..  code-block::

    cd TomBino
    mkdir build
    cd build
    cmake ..
    make
    
The cpp tests are contained in the folder `cpp_tests`, the executables can be found in `build/cpp_tests`.


To use TomBino in your code, set the compiler include path properly, and include the header files

..  code-block::

    #include <vector.h>
    #include <matrix.h>

All objects are implemented in the namespace Tombino_bla. To use them with less typing, you can set

..  code-block::
    
    namespace bla = Tombino_bla;

or even

..  code-block::
    
    using namespace Tombino_bla;

    

You can create vectors and compute with vectors like:

..  code-block:: cpp
   // creates three vectors of size 5
   Vector<double> x(5), y(5), z(5);

   // initialize the vector x with the values 0,1,2,3,4
   for (int i = 0; i < x.Size(); i++)
      x(i) = i;
   
   // initialize the vector y with the values 5,5,5,5,5
   y = 5.0;

   // compute z = x + 3*y using expression templates
   z = x+3*y;

   // print the result
   std::cout << "z = " << z << std::endl;

The class Vector is a subclass of VectorView that is a slim class that helps in the visualization of a vector.

.. code-block:: cpp
   Size() // returns the size of the vector
   Dist() // returns the distance of the vector
   View() // returns the VectorView
   Range(size_t first, size_t next) // returns a VectorView of the range [first,next)
   Slice(size_t first, size_t slice) // returns a VectorView of the sliced vector every slice elements starting from first

   operator() // access the elements of the vector
   operator= // assign a vector to another vector
   operator+ // add two vectors
   operator- // subtract two vectors
   operator* // multiply a vector by a scalar
   operator^ // multiply two vectors element by element, needs <matrix.h>



For matrices you can choose between row-major (`RowMajor`) or column-major (`ColMajor`) storage,
default is row-major.

..  code-block:: cpp

   Matrix<double,RowMajor> m1(5,3), m2(3,3);
   for (int i = 0; i < m1.Height(); i++)
     for (int j = 0; j < m1.Width(); j++)
       m1(i,j) = i+j;
   m2 = 3.7;
   Matrix product = m1 * m2;
   
You can extract a rows or a columns from a matrix:

..  code-block:: cpp

   Vector col1 = product.Col(1);

The multiplication by a scalar needs to be on the RHS:


.. code-block:: cpp

   I = dcomplex(0.0,1.0);

   Matrix<double,RowMajor> m1(5,3);
   m1 = 3.7 * m1;

   Matrix<dcomplex, RowMajor> m2(5,3);
   m2 = I * m2;


some changes ...  

   
.. toctree::
   :maxdepth: 2
   :caption: Contents:



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


