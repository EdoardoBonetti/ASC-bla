#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>
#include <string>

// include the directory when naming header files "matrix.h" and "vector.h"
#include "expression.h"
#include "lapack_interface.h"
#include "matrix.h"
#include "vector.h"

extern "C" {
#include <clapack.h>

#include "clapack.h"
}

using std::cout;
using std::endl;
using Tombino_bla::dcomplex;
using Tombino_bla::Matrix;
using Tombino_bla::Vec;
using Tombino_bla::Vector;

namespace py = pybind11;

// from ngsolve
void InitSlice(const py::slice& inds, size_t len, size_t& start, size_t& stop,
               size_t& step, size_t& n) {
  if (!inds.compute(len, &start, &stop, &step, &n))
    throw py::error_already_set();
}

////////////////////////////
// template <class T>
// class Foo {
//  public:
//   Foo(T bar) : bar_(bar) {}
//   void print() { std::cout << "Type id: " << typeid(T).name() << '\n'; }
//
//  private:
//   T bar_;
// };
////////////////////////////

// comment to be ersd

///////////////////////////
// template <typename T>
// void declare_foo(py::module& m, const std::string& typestr) {
//   using Class = Foo<T>;
//   std::string pyclass_name = std::string("Foo") + typestr;
//   py::class_<Class>(m,
//   pyclass_name.c_str()).def(py::init<T>()).def("print", &Class::print);
// }
///////////////////////////

template <typename T>
void declare_vector_class(py::module& m, const std::string& typestr) {
  using Class = Vector<T>;
  std::string pyclass_name = std::string("Vector") + typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol())
      // create a vector of given size and initialize to zero
      .def(py::init<size_t>(), py::arg("size") = 0,
           "create vector of given size and initialize to zero")

      .def("__len__", &Vector<T>::Size, "return size of vector")
      .def("__setitem__",
           [](Vector<T>& self, int i, T v) {
             if (i < 0) i += self.Size();
             self(i) = v;
           })
      .def("__getitem__",
           [](Vector<T>& self, int i) {
             if (i < 0) i += self.Size();
             return self(i);
           })
      .def("__setitem__",
           [](Vector<T>& self, py::slice inds, T val) {
             size_t start, stop, step, n;
             if (!inds.compute(self.Size(), &start, &stop, &step, &n))
               throw py::error_already_set();
             self.Range(start, stop).Slice(0, step) = val;
           })
      .def("__getitem__",
           [](Vector<T>& self, py::slice inds) {
             size_t start, stop, step, n;
             if (!inds.compute(self.Size(), &start, &stop, &step, &n))
               throw py::error_already_set();
             // return Vector<T>(self.Range(start, stop).Slice(0, step));
             return Vector(self.Range(start, stop).Slice(0, step));
           })

      .def("__add__", [](Vector<T>& self,
                         Vector<T>& other) { return Vector<T>(self + other); })

      .def("__rmul__",
           [](Vector<T>& self, T scal) { return Vector<T>(scal * self); })
      .def("__sub__", [](Vector<T>& self,
                         Vector<T>& other) { return Vector<T>(self - other); })
      .def("__rsub__", [](Vector<T>& self,
                          Vector<T>& other) { return Vector<T>(other - self); })
      .def("__mul__",
           [](Vector<T>& self, T scal) { return Vector<T>(scal * self); })
      .def("__mul__",
           [](Vector<T>& self, Vector<T>& other) { return self * other; })
      .def("__rmul__",
           [](Vector<T>& self, Vector<T>& other) { return self * other; })

      .def("__str__",
           [](const Vector<T>& self) {
             std::stringstream str;
             str << self;
             return str.str();
           })
      .def(py::pickle(
          [](Vector<T>& self) {
            // it encodes the vector as a tuple of size and bytes: it has 2
            // elements the first is the data and then size
            return py::make_tuple(

                py::bytes((char*)(void*)&self(0), self.Size() * sizeof(T)),
                self.Size());
          },
          [](py::tuple t) {
            if (t.size() != 2) throw std::runtime_error("Invalid state!");

            // Vector<T> v(size, data);
            py::bytes mem = t[0].cast<py::bytes>();
            size_t size = t[1].cast<size_t>();
            Vector<T> v(size);
            std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                        size * sizeof(T));
            return v;
          }))
      .def_buffer([](Vector<T>& v) -> py::buffer_info {
        return py::buffer_info(
            v.Data(),                           /* Pointer to buffer */
            sizeof(T),                          /* Size of one scalar */
            py::format_descriptor<T>::format(), /* Python struct-style
                                                    format descriptor */
            1,                                  /* Number of dimensions */
            {v.Size()},                         /* Buffer dimensions */
            {sizeof(T)});
      });
}

template <typename T, Tombino_bla::ORDERING ORD>
void declare_matrix_class(py::module& m, const std::string& typestr) {
  using Class = Matrix<T, ORD>;
  std::string pyclass_name = std::string("Matrix") + typestr;
  py::class_<Class>(m, pyclass_name.c_str(), py::buffer_protocol())
      .def(py::init([](size_t rows, size_t cols) {
        Matrix<T, ORD> m(rows, cols);
        if constexpr (std::is_same<T, dcomplex>::value) {
          m = dcomplex(0, 0);
        } else {
          m = 0;
        };
        return m;
      }))
      .def("__len__", &Matrix<T, ORD>::SizeRows, "return size of matrix")
      .def("__getitem__",
           [](Matrix<T, ORD>& self, std::tuple<int, int> i) {
             if (std::get<0>(i) < 0) std::get<0>(i) += self.SizeRows();
             if (std::get<1>(i) < 0) std::get<1>(i) += self.SizeCols();
             if (std::get<0>(i) < 0 ||
                 static_cast<size_t>(std::get<0>(i)) >= self.SizeRows() ||
                 std::get<1>(i) < 0 ||
                 static_cast<size_t>(std::get<1>(i)) >= self.SizeCols())
               throw py::index_error("matrix index out of range");
             return self(std::get<0>(i), std::get<1>(i));
           })
      .def("__setitem__",
           [](Matrix<T, ORD>& self, std::tuple<int, int> i, T val) {
             self(std::get<0>(i), std::get<1>(i)) = val;
           })
      .def("__getitem__",
           [](Matrix<T, ORD>& self, std::tuple<py::slice, int> inds) {
             size_t start, stop, step, n;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<1>(inds) < 0) std::get<1>(inds) += self.SizeCols();
             return Vector<T>(self.Col((size_t)std::get<1>(inds))
                                  .Range(start, stop)
                                  .Slice(0, step));
           })
      .def("__setitem__",
           [](Matrix<T, ORD>& self, std::tuple<py::slice, int> inds, T val) {
             size_t start, stop, step, n;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<1>(inds) < 0) std::get<1>(inds) += self.SizeCols();
             self.Col((size_t)std::get<1>(inds))
                 .Range(start, stop)
                 .Slice(0, step);
           })
      .def("__getitem__",
           [](Matrix<T, ORD>& self, std::tuple<int, py::slice> inds) {
             size_t start, stop, step, n;
             if (!std::get<1>(inds).compute(self.SizeCols(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<0>(inds) < 0) std::get<0>(inds) += self.SizeRows();
             return Vector<T>(self.Row((size_t)std::get<0>(inds))
                                  .Range(start, stop)
                                  .Slice(0, step));
           })
      .def("__setitem__",
           [](Matrix<T, ORD>& self, std::tuple<int, py::slice> inds, T val) {
             size_t start, stop, step, n;
             if (!std::get<1>(inds).compute(self.SizeCols(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<0>(inds) < 0) std::get<0>(inds) += self.SizeRows();
             self.Row((size_t)std::get<0>(inds))
                 .Range(start, stop)
                 .Slice(0, step) = val;
           })
      .def("__getitem__",
           [](Matrix<T, ORD>& self, std::tuple<py::slice, py::slice> inds) {
             size_t start1, stop1, step1, n1;
             size_t start2, stop2, step2, n2;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start1, &stop1,
                                            &step1, &n1))
               throw py::error_already_set();
             if (!std::get<1>(inds).compute(self.SizeCols(), &start2, &stop2,
                                            &step2, &n2))
               throw py::error_already_set();

             return Matrix<T, ORD>(self.Rows(start1, stop1)
                                       .Cols(start2, stop2)
                                       .RSlice(0, step1)
                                       .CSlice(0, step2));
           })
      .def("__setitem__",
           [](Matrix<T, ORD>& self, std::tuple<py::slice, py::slice> inds,
              T val) {
             size_t start1, stop1, step1, n1;
             size_t start2, stop2, step2, n2;

             if (!std::get<0>(inds).compute(self.SizeRows(), &start1, &stop1,
                                            &step1, &n1))
               throw py::error_already_set();
             if (!std::get<1>(inds).compute(self.SizeCols(), &start2, &stop2,
                                            &step2, &n2))
               throw py::error_already_set();
             self.Rows(start1, stop1)
                 .Cols(start2, stop2)
                 .RSlice(0, step1)
                 .CSlice(0, step2) = val;
           })
      .def("__add__",
           [](Matrix<T, ORD>& self, Matrix<T, ORD>& other) {
             return Matrix<T, ORD>(self + other);
           })
      .def("__sub__",
           [](Matrix<T, ORD>& self, Matrix<T, ORD>& other) {
             return Matrix<T, ORD>(self - other);
           })
      .def("__sub__",
           [](Matrix<T, ORD>& self, T scal) {
             Matrix<T, ORD> res(self);
             res = scal;
             return Matrix<T, ORD>(self - res);
           })
      .def("__add__",
           [](Matrix<T, ORD>& self, T scal) {
             Matrix<T, ORD> res(self);
             res = scal;
             return Matrix<T, ORD>(self + res);
           })
      .def("__rmul__", [](Matrix<T, ORD>& self,
                          T scal) { return Matrix<T, ORD>(scal * self); })
      .def("__mul__", [](Matrix<T, ORD>& self,
                         T scal) { return Matrix<T, ORD>(scal * self); })
      .def("__mul__",
           [](Matrix<T, ORD>& self, Matrix<T, ORD>& other) {
             return Matrix<T, ORD>(self * other);
           })
      .def("__rmul__",
           [](Matrix<T, ORD>& self, Matrix<T, ORD>& other) {
             return Matrix<T, ORD>(self * other);
           })
      .def("__mul__", [](Matrix<T, ORD>& self,
                         Vector<T>& other) { return Vector<T>(self * other); })
      .def("__rmul__", [](Matrix<T, ORD>& self,
                          Vector<T>& other) { return Vector<T>(self * other); })

      .def("__str__",
           [](const Matrix<T, ORD>& self) {
             std::stringstream str;
             str << self;
             return str.str();
           })
      .def_property_readonly("shape",
                             [](const Matrix<T, ORD>& self) {
                               return std::tuple(self.SizeRows(),
                                                 self.SizeCols());
                             })
      .def_property_readonly("T",
                             [](const Matrix<T, ORD>& self) {
                               return Matrix<T, ORD>(Transpose(self));
                             })
      .def_property_readonly("inv",
                             [](const Matrix<T, ORD>& self) {
                               return Matrix<T, ORD>(Inverse(self));
                             })

      .def(py::pickle(
          [](Matrix<T, ORD>& self) {
            return py::make_tuple(
                self.SizeRows(), self.SizeCols(),
                py::bytes((char*)(void*)&self(0, 0),
                          self.SizeRows() * self.SizeCols() * sizeof(T)));
          },
          [](py::tuple t) {
            if (t.size() != 2) throw std::runtime_error("should be a 2-tuple!");

            Matrix<T, ORD> v(t[0].cast<size_t>(), t[1].cast<size_t>());
            py::bytes mem = t[2].cast<py::bytes>();
            std::memcpy(&v(0, 0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                        v.SizeCols() * v.SizeRows() * sizeof(T));
            return v;
          }))

      .def_buffer([](Matrix<T, ORD>& v) -> py::buffer_info {
        return py::buffer_info(
            v.Data(),                           /* Pointer to buffer */
            sizeof(T),                          /* Size of one scalar */
            py::format_descriptor<T>::format(), /* Python struct-style
                                                        format descriptor */
            2,                                  /* Number of dimensions */
            {v.SizeRows(), v.SizeCols()},       /* Buffer dimensions */
            {sizeof(T) * v.SizeCols(), /* Strides (in bytes) for each index */
             sizeof(T)});
      });
}

template <int SIZE, typename T>
void declare_vec_class(py::module& m, const std::string& typestr) {
  using Class = Vec<SIZE, T>;
  std::string pyclass_name = std::string("Vec") + typestr;
  py::class_<Class>(m, pyclass_name.c_str())
      .def(py::init<>(), "create Vec")
      .def("__len__", &Vec<SIZE, T>::Size, "return size of vector")
      .def("__setitem__",
           [](Vec<SIZE, T>& self, int i, T v) {
             if (i < 0) i += self.Size();
             if (i < 0 || static_cast<size_t>(i) >= self.Size())
               throw py::index_error("vector index out of range");
             self(i) = v;
           })
      .def("__getitem__",
           [](Vec<SIZE, T>& self, int i) {
             if (i < 0) i += self.Size();
             if (i < 0 || static_cast<size_t>(i) >= self.Size())
               throw py::index_error("vector index out of range");
             return self(i);
           })
      .def("__str__",
           [](const Vec<SIZE, T>& self) {
             std::stringstream str;
             str << self;
             return str.str();
           })

      //      .def("__add__", [](Vector<T>& self, Vector<T>& other)
      //           { return Vector<T>(self + other); })
      //
      //      .def("__rmul__",
      //           [](Vector<T>& self, T scal) { return Vector<T>(scal * self);
      //           })
      //
      //      .def("__mul__",
      //           [](Vector<T>& self, T scal) { return Vector<T>(scal * self);
      //           })
      //      .def("__mul__",
      //           [](Vector<T>& self, Vector<T>& other) { return self * other;
      //           })
      //      .def("__rmul__",
      //           [](Vector<T>& self, Vector<T>& other) { return self * other;
      //           })
      // .def("__str__",
      //     [](const Vec<SIZE, T>& self)
      //     {
      //       std::stringstream str;
      //       str << self;
      //       return str.str();
      //     })

      .def(py::pickle(
          [](Vec<SIZE, T>& self) {
            return py::make_tuple(
                self.Size(),
                py::bytes((char*)(void*)&self(0), self.Size() * sizeof(T)));
          },
          [](py::tuple t) {
            if (t.size() != 2) throw std::runtime_error("should be a 2-tuple!");
            Vec<SIZE, T> v(t[0].cast<size_t>());
            py::bytes mem = t[1].cast<py::bytes>();
            std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                        v.Size() * sizeof(T));
            return v;
          }));
}

PYBIND11_MODULE(bla, m) {
  m.doc() = "Basic linear algebra module";  // optional module docstring
  declare_vector_class<int>(m, "Int");
  declare_vector_class<double>(m, "");
  // declare_vector_class<dcomplex>(m, "Complex");

  declare_matrix_class<double, Tombino_bla::ORDERING::RowMajor>(m, "");
  declare_matrix_class<dcomplex, Tombino_bla::ORDERING::RowMajor>(m, "Complex");

  declare_vec_class<2, double>(m, "2");
  declare_vec_class<3, double>(m, "3");

  // define the inner product
  m.def("InnerProduct",
        [](Matrix<double, Tombino_bla::ORDERING::RowMajor>& self,
           Matrix<double, Tombino_bla::ORDERING::RowMajor>& other) {
          return InnerProduct(self, other);
        });

  m.def("InnerProduct", [](Vector<double>& self, Vector<double>& other) {
    return InnerProduct(self, other);
  });
  m.def("Inverse", [](Matrix<double, Tombino_bla::ORDERING::RowMajor>& self) {
    return Inverse(self);
  });
  m.def("Transpose", [](Matrix<double, Tombino_bla::ORDERING::RowMajor>& self) {
    return Matrix<double, Tombino_bla::ORDERING::RowMajor>(Transpose(self));
  });

  // declare_matrix_class<double, ColMajor>(m, "ColMajor");
  // declare_matrix_class<dcomplex, ColMajor>(m, "ComplexColMajor");

  // and 4 classes : LapackLU, LapackQR,LapackEVP,  LapackSVD.

  // py::class_<LapackLU<RowMajor>>(m, "LapackLU")
  //     .def(py::init<Matrix<double, RowMajor>>(), py::arg("A"), "create LU
  //     decomposition of given matrix") .def("__str__",
  //          [](const LapackLU<RowMajor>& self) {
  //            std::stringstream str;
  //            str << self;
  //            return str.str();
  //          })
  //     .def("L", [](const LapackLU<RowMajor>& self) { self.LFactor(); });

  // we wrap the lapack mult matrix matrix

  // now we wrap again the vector class but using a template
}
