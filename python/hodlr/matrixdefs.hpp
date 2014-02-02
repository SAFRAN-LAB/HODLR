#ifndef __MATRIXDEFS_HPP__
#define __MATRIXDEFS_HPP__

#include <Python.h>
#include "HODLR_Matrix.hpp"

class Python_Matrix : public HODLR_Matrix {

public:

    Python_Matrix (PyObject * function) : function_(function) {
        status_ = 0;
        Py_XINCREF(function_);
    };

    ~Python_Matrix () {
        Py_XDECREF(function_);
    };

    double get_Matrix_Entry (const unsigned i, const unsigned j) {
        double value;
        PyObject *arglist, *result;

        // Set up the arguments for the Python function.
        arglist = Py_BuildValue("ii", i, j);
        if (arglist == NULL) {
            Py_XDECREF(arglist);
            status_ = 1;
            return 0.0;
        }

        // Call the Python function.
        result = PyObject_CallObject(function_, arglist);
        Py_DECREF(arglist);

        // Parse the result or fail.
        if (result == NULL || !PyFloat_Check(result)) {
            Py_XDECREF(result);
            status_ = 2;
            return 0.0;
        }
        value = PyFloat_AsDouble(result);

        // Clean up and return.
        Py_DECREF(result);
        return value;
    };

private:

    unsigned int status_;
    PyObject * function_;

};

#endif /* defined(__MATRIXDEFS_HPP__) */
