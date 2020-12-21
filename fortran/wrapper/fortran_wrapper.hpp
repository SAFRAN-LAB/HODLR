#ifndef __fortran_wrapper__
#define __fortran_wrapper__

#include "Kernel.hpp"
#include "HODLR/LowRank.hpp"
#include "HODLR/HODLR.hpp"

// Extern C needs to be used to prevent name mangling
extern "C"
{
    void initialize_kernel_object_c(int N, Kernel** kernel);
    void get_matrix_c(Kernel** kernel, int row_start, int col_start, int row_end, int col_end, double* matrix);
    void initialize_lowrank_c(Kernel** kernel, char* factorization_type, LowRank** lowrank);
    void get_rank_c(LowRank** lowrank, double eps, int& rank);
    void get_lr_c(LowRank** lowrank, double* l, double* r);
    void initialize_hodlr_c(int N, int M, double eps, HODLR** tree);
    void assemble_c(HODLR** tree, Kernel** kernel, char* lowrank_method, bool is_sym, bool is_pd);
    void matmat_product_c(HODLR** tree, double* x, double* b);
    void factorize_c(HODLR** tree);
    void solve_c(HODLR** tree, double* b, double* x);
    void logdeterminant_c(HODLR** tree, double &log_det);
    void symm_factor_transpose_prod_c(HODLR** tree, double* x, double* b);
    void symm_factor_prod_c(HODLR** tree, double* x, double* b);
    void get_symm_factor_c(HODLR** tree, double* W_matrix);
}

#endif /* (defined(__fortran_wrapper__) */
