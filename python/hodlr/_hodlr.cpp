#include <Python.h>
#include <structmember.h>
#include <numpy/arrayobject.h>
#include <Eigen/Dense>

#include "HODLR_Tree.hpp"
#include "matrixdefs.hpp"

// Silence some warnings.
#pragma GCC diagnostic ignored "-Wwrite-strings"
#pragma GCC diagnostic ignored "-Wdelete-non-virtual-dtor"

using namespace Eigen;

// Parse NumPy objects.
#define PARSE_ARRAY(o) (PyArrayObject*) PyArray_FROM_OTF(o, NPY_DOUBLE, \
        NPY_IN_ARRAY)

// The Python object type.
typedef struct {
    PyObject_HEAD

    unsigned int dim;
    Gaussian_Matrix * matrix;
    HODLR_Tree<Gaussian_Matrix> * solver;
    PyArrayObject * time_array;

} _hodlr;

// Avoid name mangling.
extern "C" {
    static void _hodlr_dealloc(_hodlr *self);
    static PyObject *_hodlr_new(PyTypeObject *type, PyObject *args, PyObject *kwds);
    static int _hodlr_init(_hodlr *self, PyObject *args, PyObject *kwds);
    void init_hodlr(void);

    static PyObject* _hodlr_logdet (_hodlr *self);
    static PyObject* _hodlr_solve (_hodlr *self, PyObject *args);
    static PyObject* _hodlr_matrix_product (_hodlr *self, PyObject *args);
}

// General Python setup and cleanup.
static void _hodlr_dealloc(_hodlr *self)
{
    if (self->solver != NULL) delete self->solver;
    if (self->matrix != NULL) delete self->matrix;
    Py_XDECREF(self->time_array);
    self->ob_type->tp_free((PyObject*)self);
}

static PyObject *_hodlr_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    _hodlr* self;
    self = (_hodlr*)type->tp_alloc(type, 0);
    self->matrix = NULL;
    self->solver = NULL;
    self->time_array = NULL;
    return (PyObject*)self;
}

static int _hodlr_init(_hodlr *self, PyObject *args, PyObject *kwds)
{
    double tol = 1e-14, amp, var;
    unsigned int nLeaf = 50;
    PyObject * time_obj = NULL, * diag_obj = NULL;
    static char *kws[] = {"amp", "var", "time", "diag", "nleaf", "tol", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddOO|id", &(kws[0]), &amp,
                                     &var, &time_obj, &diag_obj, &nLeaf,
                                     &tol))
        return -1;

    // Set up the matrix.
    self->time_array = PARSE_ARRAY(time_obj);
    self->dim = PyArray_DIM(self->time_array, 0);
    double *time = (double*)PyArray_DATA(self->time_array);
    if (self->time_array == NULL) {
        return -2;
    }
    self->matrix = new Gaussian_Matrix (amp, var, time);

    // Parse the diagonal array.
    PyArrayObject *diag_array = PARSE_ARRAY(diag_obj);
    if (diag_array == NULL || self->dim != PyArray_DIM(diag_array, 0)) {
        Py_DECREF(self->time_array);
        Py_XDECREF(diag_array);
        return -3;
    }
    double *diag = (double*)PyArray_DATA(diag_array);
    VectorXd dv = Map<VectorXd>(diag, self->dim);

    // Initialize the solver.
    self->solver = new HODLR_Tree<Gaussian_Matrix> (self->matrix, self->dim, nLeaf);
    self->solver->assemble_Matrix(dv, tol);
    Py_DECREF(diag_array);

    // Factorize the matrix.
    self->solver->compute_Factor();

    return 0;
}

/* static int _hodlr_init(_hodlr *self, PyObject *args, PyObject *kwds) */
/* { */
/*     double tol = 1e-14; */
/*     unsigned int nLeaf = 50; */
/*     PyObject * matrix_obj = NULL, * diag_obj = NULL; */
/*     static char *kws[] = {"matrix", "diag", "nleaf", "tol", NULL}; */
/*     if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|id", &(kws[0]), */
/*                                      &matrix_obj, &diag_obj, &nLeaf, &tol)) */
/*         return -1; */

/*     // Check the matrix function and set up the matrix object. */
/*     if (!PyCallable_Check (matrix_obj)) { */
/*         PyErr_SetString(PyExc_TypeError, "The matrix object must be callable"); */
/*         return -2; */
/*     } */
/*     self->matrix = new Python_Matrix (matrix_obj); */

/*     // Parse the diagonal array. */
/*     PyArrayObject *diag_array = PARSE_ARRAY(diag_obj); */
/*     if (diag_array == NULL) { */
/*         Py_XDECREF(diag_array); */
/*         return -3; */
/*     } */
/*     self->dim = PyArray_DIM(diag_array, 0); */
/*     double *diag = (double*)PyArray_DATA(diag_array); */
/*     VectorXd dv = Map<VectorXd>(diag, self->dim); */

/*     // Initialize the solver. */
/*     self->solver = new HODLR_Tree<Python_Matrix> (self->matrix, self->dim, nLeaf); */
/*     self->solver->assemble_Matrix(dv, tol); */
/*     Py_DECREF(diag_obj); */

/*     // Factorize the matrix. */
/*     self->solver->compute_Factor(); */

/*     return 0; */
/* } */

static PyObject *_hodlr_logdet (_hodlr *self)
{
    double logdet;
    self->solver->compute_Determinant(logdet);
    return Py_BuildValue("d", logdet);
}

static PyObject *_hodlr_solve (_hodlr *self, PyObject *args)
{
    PyObject *b_obj;
    if (!PyArg_ParseTuple(args, "O", &b_obj)) return NULL;

    // Parse the RHS array.
    PyArrayObject *b_array = PARSE_ARRAY(b_obj);
    if (b_array == NULL) {
        Py_XDECREF(b_array);
        return NULL;
    }

    // Check the dimensions.
    unsigned int ndim = PyArray_NDIM(b_array);
    if (ndim < 1 || ndim > 2 || self->dim != PyArray_DIM(b_array, 0)) {
        PyErr_SetString(PyExc_ValueError, "Dimension mismatch");
        Py_DECREF(b_array);
        return NULL;
    }

    // How many systems are we solving?
    unsigned int nrhs = 1, oned = 1;
    if (ndim > 1) {
        oned = 0;
        nrhs = PyArray_DIM(b_array, 1);
    }

    // Access the data.
    double *b = (double*)PyArray_DATA(b_array);
    MatrixXd bv = Map<MatrixXd>(b, self->dim, nrhs);

    // Solve the system.
    MatrixXd xv(self->dim, nrhs);
    self->solver->solve(bv, xv);
    Py_DECREF(b_array);

    // Build the output array.
    PyArrayObject *x_array;
    if (oned) {
        npy_intp dim[1] = {self->dim};
        x_array = (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    } else {
        npy_intp dim[2] = {self->dim, nrhs};
        x_array = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    }
    if (x_array == NULL) {
        Py_DECREF(b_array);
        Py_XDECREF(x_array);
        return NULL;
    }
    int i, j;
    double *x = (double*)PyArray_DATA(x_array);
    for (i = 0; i < self->dim; ++i)
        for (j = 0; j < nrhs; ++j)
            x[i*nrhs+j] = xv(i, j);

    // Build the return value.
    PyObject *ret = Py_BuildValue("O", x_array);
    Py_DECREF(x_array);

    if (ret == NULL) {
        Py_XDECREF(ret);
        return NULL;
    }

    return ret;
}

static PyObject *_hodlr_matrix_product (_hodlr *self, PyObject *args)
{
    PyObject *b_obj;
    if (!PyArg_ParseTuple(args, "O", &b_obj)) return NULL;

    // Parse the matrix.
    PyArrayObject *b_array = PARSE_ARRAY(b_obj);
    if (b_array == NULL) {
        Py_XDECREF(b_array);
        return NULL;
    }

    // Check the dimensions.
    unsigned int ndim = PyArray_NDIM(b_array);
    if (ndim < 1 || ndim > 2 || self->dim != PyArray_DIM(b_array, 0)) {
        PyErr_SetString(PyExc_ValueError, "Dimension mismatch");
        Py_DECREF(b_array);
        return NULL;
    }

    // How many systems are we solving?
    unsigned int nrhs = 1, oned = 1;
    if (ndim > 1) {
        oned = 0;
        nrhs = PyArray_DIM(b_array, 1);
    }

    // Access the data.
    double *b = (double*)PyArray_DATA(b_array);
    MatrixXd bv = Map<MatrixXd>(b, self->dim, nrhs);

    // Solve the system.
    MatrixXd xv(self->dim, nrhs);
    self->solver->matMatProduct(bv, xv);
    Py_DECREF(b_array);

    // Build the output array.
    PyArrayObject *x_array;
    if (oned) {
        npy_intp dim[1] = {self->dim};
        x_array = (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);
    } else {
        npy_intp dim[2] = {self->dim, nrhs};
        x_array = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    }
    if (x_array == NULL) {
        Py_DECREF(b_array);
        Py_XDECREF(x_array);
        return NULL;
    }
    int i, j;
    double *x = (double*)PyArray_DATA(x_array);
    for (i = 0; i < self->dim; ++i)
        for (j = 0; j < nrhs; ++j)
            x[i*nrhs+j] = xv(i, j);

    // Build the return value.
    PyObject *ret = Py_BuildValue("O", x_array);
    Py_DECREF(x_array);

    if (ret == NULL) {
        Py_XDECREF(ret);
        return NULL;
    }

    return ret;
}

static PyMethodDef _hodlr_methods[] = {
    {"logdet",
     (PyCFunction)_hodlr_logdet,
     METH_NOARGS,
     "Get the log-determinant of the system."},
    {"solve",
     (PyCFunction)_hodlr_solve,
     METH_VARARGS,
     "Solve the system for a given RHS."},
    {"matrix_product",
     (PyCFunction)_hodlr_matrix_product,
     METH_VARARGS,
     "Perform a matrix multiplication."},
    {NULL}  /* Sentinel */
};

static PyMemberDef _hodlr_members[] = {{NULL}};

static char _hodlr_doc[] = "This is the ``_hodlr`` object. There is some black magic.";
static PyTypeObject _hodlr_type = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_hodlr.HODLR",            /*tp_name*/
    sizeof(_hodlr),            /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)_hodlr_dealloc,/*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    _hodlr_doc,                /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    _hodlr_methods,            /* tp_methods */
    _hodlr_members,            /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    (initproc)_hodlr_init,     /* tp_init */
    0,                         /* tp_alloc */
    _hodlr_new,                /* tp_new */
};


//
// Initialize the module.
//

static char module_doc[] = "HODLR Solver";
static PyMethodDef module_methods[] = {{NULL}};
void init_hodlr(void)
{
    PyObject *m;

    if (PyType_Ready(&_hodlr_type) < 0)
        return;

    m = Py_InitModule3("_hodlr", module_methods, module_doc);
    if (m == NULL)
        return;

    Py_INCREF(&_hodlr_type);
    PyModule_AddObject(m, "HODLR", (PyObject *)&_hodlr_type);

    import_array();
}
