#include "HODLR_Matrix.hpp"
#include "Matrix_Factorizer.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;

class PyMatrix : public HODLR_Matrix
{
public:
    /* Inherit the constructors */
    using HODLR_Matrix::HODLR_Matrix;

    /* Trampoline (need one for each virtual function) */
    dtype getMatrixEntry(int i, int j) override 
    {
        PYBIND11_OVERLOAD_PURE(dtype,          /* Return type */
                               HODLR_Matrix,   /* Parent class */
                               getMatrixEntry, /* Name of function in C++ (must match Python name) */
                               i, j            /* Arguments */
                              );
    }
};

PYBIND11_MODULE(pyhodlrlib, m) 
{
    m.doc() = "This is the Python Wrapper to HODLRlib";

    py::class_<HODLR_Matrix, PyMatrix> hodlr_matrix(m, "HODLR_Matrix");
    hodlr_matrix
        .def(py::init<int>())
        .def("getMatrixEntry", &HODLR_Matrix::getMatrixEntry)
        .def("getMatrix", &HODLR_Matrix::getMatrix);

    py::class_<Matrix_Factorizer> matrix_factorizer(m, "Matrix_Factorizer");
    matrix_factorizer
        .def(py::init<HODLR_Matrix*, std::string>())
        .def("factorize", &Matrix_Factorizer::factorize)
        .def("getL", &Matrix_Factorizer::getL)
        .def("getR", &Matrix_Factorizer::getR);
}
