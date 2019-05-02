#ifndef __fortran_wrapper__
#define __fortran_wrapper__

#include "Kernel.hpp"
#include "Matrix_Factorizer.hpp"
#include "HODLR_Tree.hpp"

// Extern C needs to be used to prevent name mangling
extern "C"
{
    void initialize_kernel_object_c(Kernel** kernel, int N);
    void get_matrix_c(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end);
    void initialize_matrix_factorizer_c(Matrix_Factorizer** factorizer, Kernel** kernel, char* factorization_type);
    void get_factorization_c(Matrix_Factorizer** factorizer, double* l, double* r, double eps);
    void initialize_hodlr_tree_c(HODLR_Tree** tree, int n_levels, double eps, Matrix_Factorizer** factorizer);
    void assemble_tree_c(HODLR_Tree** tree, bool is_sym, bool is_pd);
    void matmat_product_c(HODLR_Tree** tree, double* x, double* b);
    void factorize_c(HODLR_Tree** tree);
    void solve_c(HODLR_Tree** tree, double* b, double* x);
    void logdeterminant_c(HODLR_Tree** tree, double &log_det);
}

#endif /* (defined(__fortran_wrapper__) */
