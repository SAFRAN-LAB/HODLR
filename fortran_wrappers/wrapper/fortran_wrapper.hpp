#ifndef __fortran_wrapper__
#define __fortran_wrapper__

#include "Kernel.hpp"
#include "Matrix_Factorizer.hpp"
// #include "HODLR_Tree.hpp"

// Extern C needs to be used to prevent name mangling
extern "C"
{
    void initialize_kernel_object(Kernel** kernel, int N, int dim);
    void get_matrix(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end);
    void initialize_matrix_factorizer(Matrix_Factorizer** factorizer, Kernel** kernel, char* factorization_type);
    void get_factorization(Matrix_Factorizer** factorizer, double* l, double* r, double eps);
    // void initialize_hodlr_tree(HODLR_Tree** tree, int n_levels, double eps, Matrix_Factorizer** factorizer);
}

#endif /* (defined(__fortran_wrapper__) */
