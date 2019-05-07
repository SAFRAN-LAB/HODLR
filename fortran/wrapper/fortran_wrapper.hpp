#ifndef __fortran_wrapper__
#define __fortran_wrapper__

#include "Kernel.hpp"
#include "Matrix_Factorizer.hpp"
#include "HODLR_Tree.hpp"

// Extern C needs to be used to prevent name mangling
extern "C"
{
    void initialize_kernel_object_c(int N, Kernel** kernel);
    void get_matrix_c(Kernel** kernel, int row_start, int col_start, int row_end, int col_end, double* matrix);
    void initialize_matrix_factorizer_c(Kernel** kernel, char* factorization_type, Matrix_Factorizer** factorizer);
    void get_factorization_c(Matrix_Factorizer** factorizer, double eps, double* l, double* r, int& rank);
    void initialize_hodlr_tree_c(int n_levels, double eps, Matrix_Factorizer** factorizer, HODLR_Tree** tree);
    void assemble_tree_c(HODLR_Tree** tree, bool is_sym, bool is_pd);
    void matmat_product_c(HODLR_Tree** tree, double* x, double* b);
    void factorize_c(HODLR_Tree** tree);
    void solve_c(HODLR_Tree** tree, double* b, double* x);
    void logdeterminant_c(HODLR_Tree** tree, double &log_det);
    void symm_factor_transpose_prod_c(HODLR_Tree** tree, double* x, double* b);
    void symm_factor_prod_c(HODLR_Tree** tree, double* x, double* b);
    void get_symm_factor_c(HODLR_Tree** tree, double* W_matrix);
}

#endif /* (defined(__fortran_wrapper__) */
