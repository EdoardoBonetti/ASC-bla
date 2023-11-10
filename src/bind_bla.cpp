#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/pybind11.h>

#include <sstream>

#include "matrix.h"
#include "vector.h"
// #include "matrix.h"

#include "lapack_interface.h"

using namespace Tombino_bla;
namespace py = pybind11;

extern "C" {
#include <clapack.h>

#include "clapack.h"
}

PYBIND11_MODULE(bla, m) {
  m.doc() = "Basic linear algebra module";  // optional module docstring

  py::class_<Vector<double>>(m, "Vector")
      .def(py::init<size_t>(), py::arg("size"), "create vector of given size")
      .def("__len__", &Vector<double>::Size, "return size of vector")
      .def("__setitem__",
           [](Vector<double>& self, int i, double v) {
             if (i < 0) i += self.Size();
             if (i < 0 || i >= self.Size())
               throw py::index_error("vector index out of range");
             self(i) = v;
           })
      .def("__getitem__", [](Vector<double>& self, int i) { return self(i); })

      .def("__setitem__",
           [](Vector<double>& self, py::slice inds, double val) {
             size_t start, stop, step, n;
             if (!inds.compute(self.Size(), &start, &stop, &step, &n))
               throw py::error_already_set();
             self.Range(start, stop).Slice(0, step) = val;
           })
      .def("__getitem__",
           [](Vector<double>& self, py::slice inds) {
             size_t start, stop, step, n;
             if (!inds.compute(self.Size(), &start, &stop, &step, &n))
               throw py::error_already_set();
             return self.Range(start, stop).Slice(0, step);
           })

      .def("__add__",
           [](Vector<double>& self, Vector<double>& other) {
             return Vector<double>(self + other);
           })

      .def("__rmul__", [](Vector<double>& self,
                          double scal) { return Vector<double>(scal * self); })

      .def("__mul__", [](Vector<double>& self,
                         double scal) { return Vector<double>(scal * self); })

      // Inner Product between two vectors
      .def("__mul__", [](Vector<double>& self,
                         Vector<double>& other) { return self * other; })
      .def("__rmul__", [](Vector<double>& self,
                          Vector<double>& other) { return self * other; })

      .def("__str__",
           [](const Vector<double>& self) {
             std::stringstream str;
             str << self;
             return str.str();
           })

      .def(py::pickle(
          [](Vector<double>& self) {  // __getstate__
            /* return a tuple that fully encodes the state of the object */
            return py::make_tuple(self.Size(),
                                  py::bytes((char*)(void*)&self(0),
                                            self.Size() * sizeof(double)));
          },
          [](py::tuple t) {  // __setstate__
            if (t.size() != 2) throw std::runtime_error("should be a 2-tuple!");

            Vector<double> v(t[0].cast<size_t>());
            py::bytes mem = t[1].cast<py::bytes>();
            std::memcpy(&v(0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                        v.Size() * sizeof(double));
            return v;
          }))

      ;

  // from here my Matrix implementation
  py::class_<Matrix<double, RowMajor>>(m, "Matrix", py::buffer_protocol())
      //.def(py::init<size_t, size_t>(), py::arg("rows"), py::arg("cols"),
      //     "create matrix : rows x cols")
      //  .def("__getitem__",
      //       [](Matrix<double, RowMajor> self, py::tuple i) {
      //         return self(i[0].cast<size_t>(), i[1].cast<size_t>());
      //       })
      .def(py::init([](size_t rows, size_t cols) {
        Matrix<double, RowMajor> m(rows, cols);
        m = 0;
        return m;
      }))

      // setter and getter for the single elements in the matrix
      .def("__getitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<int, int> i) {
             return self(std::get<0>(i), std::get<1>(i));
           })
      .def("__setitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<int, int> i,
              double val) { self(std::get<0>(i), std::get<1>(i)) = val; })

      // getter for the rows
      .def("__getitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<int, py::slice> inds) {
             size_t start, stop, step, n;
             if (!std::get<1>(inds).compute(self.SizeCols(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<0>(inds) < 0) std::get<0>(inds) += self.SizeRows();
             // now cast hardcast to size_t
             return self.Row((size_t)std::get<0>(inds))
                 .Range(start, stop)
                 .Slice(0, step);
           })

      // setter for the rows
      .def("__setitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<int, py::slice> inds,
              double val) {
             size_t start, stop, step, n;
             if (!std::get<1>(inds).compute(self.SizeCols(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<0>(inds) < 0) std::get<0>(inds) += self.SizeRows();
             // now cast hardcast to size_t
             self.Row((size_t)std::get<0>(inds))
                 .Range(start, stop)
                 .Slice(0, step) = val;
           })

      // getter for the cols
      .def("__getitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<py::slice, int> inds) {
             size_t start, stop, step, n;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<1>(inds) < 0) std::get<1>(inds) += self.SizeCols();
             // now cast hardcast to size_t
             return self.Col((size_t)std::get<1>(inds))
                 .Range(start, stop)
                 .Slice(0, step);
           })

      // setter for the cols
      .def("__setitem__",
           [](Matrix<double, RowMajor>& self, std::tuple<py::slice, int> inds,
              double val) {
             size_t start, stop, step, n;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start, &stop,
                                            &step, &n))
               throw py::error_already_set();
             if (std::get<1>(inds) < 0) std::get<1>(inds) += self.SizeCols();
             // now cast hardcast to size_t
             self.Col((size_t)std::get<1>(inds))
                 .Range(start, stop)
                 .Slice(0, step);
           })

      // getter for the submatrix
      .def("__getitem__",
           [](Matrix<double, RowMajor>& self,
              std::tuple<py::slice, py::slice> inds) {
             size_t start1, stop1, step1, n1;
             size_t start2, stop2, step2, n2;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start1, &stop1,
                                            &step1, &n1))
               throw py::error_already_set();
             if (!std::get<1>(inds).compute(self.SizeCols(), &start2, &stop2,
                                            &step2, &n2))
               throw py::error_already_set();
             return Matrix<double, RowMajor>(
                 self.Cols(start1, stop1).Rows(start2, stop2));
           })

      .def("__setitem__",
           [](Matrix<double, RowMajor>& self,
              std::tuple<py::slice, py::slice> inds, double val) {
             size_t start1, stop1, step1, n1;
             size_t start2, stop2, step2, n2;
             if (!std::get<0>(inds).compute(self.SizeRows(), &start1, &stop1,
                                            &step1, &n1))
               throw py::error_already_set();
             if (!std::get<1>(inds).compute(self.SizeCols(), &start2, &stop2,
                                            &step2, &n2))
               throw py::error_already_set();
             self.Cols(start1, stop1).Rows(start2, stop2) = val;
           })

      .def("__add__",
           [](Matrix<double, RowMajor>& self, Matrix<double, RowMajor>& other) {
             return Matrix<double, RowMajor>(self + other);
           })

      .def("__sub__",
           [](Matrix<double, RowMajor>& self, Matrix<double, RowMajor>& other) {
             return Matrix<double, RowMajor>(self - other);
           })
      .def("__sub__",
           [](Matrix<double, RowMajor>& self, double scal) {
             Matrix<double, RowMajor> res(self);
             res = scal;
             return Matrix<double, RowMajor>(self - res);
           })

      // addition by scalar matrix
      .def("__add__",
           [](Matrix<double, RowMajor>& self, double scal) {
             Matrix<double, RowMajor> res(self);
             res = scal;
             std::cout << "res = " << res << std::endl;
             return Matrix<double, RowMajor>(self + res);
           })

      // mult by scalar matrix
      .def("__rmul__",
           [](Matrix<double, RowMajor>& self, double scal) {
             return Matrix<double, RowMajor>(scal * self);
           })
      .def("__mul__",
           [](Matrix<double, RowMajor>& self, double scal) {
             return Matrix<double, RowMajor>(scal * self);
           })

      // mult matrix matrix
      .def("__mul__",
           [](Matrix<double, RowMajor>& self, Matrix<double, RowMajor>& other) {
             return Matrix<double, RowMajor>(self * other);
           })
      .def("__rmul__",
           [](Matrix<double, RowMajor>& self, Matrix<double, RowMajor>& other) {
             return Matrix<double, RowMajor>(self * other);
           })

      // mult matrix vector
      .def("__mul__",
           [](Matrix<double, RowMajor>& self, Vector<double>& other) {
             return Vector<double>(self * other);
           })
      .def("__rmul__",
           [](Matrix<double, RowMajor>& self, Vector<double>& other) {
             return Vector<double>(self * other);
           })

      .def("__str__",
           [](const Matrix<double, RowMajor>& self) {
             std::stringstream str;
             str << self;
             return str.str();
           })

      // shape of the matrix
      .def_property_readonly("shape",
                             [](const Matrix<double, RowMajor>& self) {
                               return std::tuple(self.SizeRows(),
                                                 self.SizeCols());
                             })

      // transpose
      .def_property_readonly("T",
                             [](const Matrix<double, RowMajor>& self) {
                               return Matrix<double, RowMajor>(Transpose(self));
                             })

      // inverse
      .def_property_readonly("inv",
                             [](const Matrix<double, RowMajor>& self) {
                               return Matrix<double, RowMajor>(Inverse(self));
                             })

      .def(py::pickle(
          [](Matrix<double, RowMajor>& self) {  // __getstate__
            /* return a tuple that fully encodes the state of the object */
            return py::make_tuple(
                self.SizeRows(), self.SizeCols(),
                py::bytes((char*)(void*)&self(0, 0),
                          self.SizeRows() * self.SizeCols() * sizeof(double)));
          },
          [](py::tuple t) {  // __setstate__
            if (t.size() != 2) throw std::runtime_error("should be a 2-tuple!");

            Matrix<double, RowMajor> v(t[0].cast<size_t>(),
                                       t[1].cast<size_t>());
            py::bytes mem = t[2].cast<py::bytes>();
            std::memcpy(&v(0, 0), PYBIND11_BYTES_AS_STRING(mem.ptr()),
                        v.SizeCols() * v.SizeRows() * sizeof(double));
            return v;
          }))

      .def_buffer([](Matrix<double, RowMajor>& v) -> py::buffer_info {
        return py::buffer_info(
            v.Data(),                                /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style
                                                        format descriptor */
            2,                                       /* Number of dimensions */
            {v.SizeRows(), v.SizeCols()},            /* Buffer dimensions */
            {sizeof(double) *
                 v.SizeCols(), /* Strides (in bytes) for each index */
             sizeof(double)});
      });

  // and 4 classes : LapackLU, LapackQR,LapackEVP,  LapackSVD.

  // py::class_<LapackLU<RowMajor>>(m, "LapackLU")
  //     .def(py::init<Matrix<double, RowMajor>>(), py::arg("A"), "create LU decomposition of given matrix")
  //     .def("__str__",
  //          [](const LapackLU<RowMajor>& self) {
  //            std::stringstream str;
  //            str << self;
  //            return str.str();
  //          })
  //     .def("L", [](const LapackLU<RowMajor>& self) { self.LFactor(); });

  // we wrap the lapack mult matrix matrix
}

