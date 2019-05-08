#include "fortran_wrapper.hpp"

// This is a helper function:
void removeEndSpaces(char *str) 
{ 
    int i;
    for (i = 0; str[i]; i++) 
        if(str[i] == ' ')
            break;

    str[i] = '\0'; 
} 

void initialize_kernel_object_c(int N, Kernel** kernel)
{
    (*kernel) = new Kernel(N);
    return;
}

void get_matrix_c(Kernel** kernel, int row_start, int col_start, int row_end, int col_end, double* matrix)
{
    Mat temp = (*kernel)->getMatrix(row_start, col_start, row_end, col_end);
    // Counter variables:
    int i, j;
    for(i = 0; i < temp.cols(); i++)
        for(j = 0; j < temp.rows(); j++)
            matrix[temp.rows() * i + j] = temp(j, i);

    return;
}

void initialize_lowrank_c(Kernel** kernel, char* factorization_type, LowRank** lowrank)
{   
    removeEndSpaces(factorization_type);
    (*lowrank) = new LowRank((*kernel), factorization_type);
    return;
}

void get_rank_c(LowRank** lowrank, double eps, int& rank)
{
    (*lowrank)->factorize(eps);
    rank = ((*lowrank)->L).cols();
}

void get_lr_c(LowRank** lowrank, double* l, double* r)
{
    // Counter variables:
    int i, j;
    for(i = 0; i < ((*lowrank)->L).cols(); i++)
        for(j = 0; j < ((*lowrank)->L).rows(); j++)
            l[((*lowrank)->L).rows() * i + j] = ((*lowrank)->L)(j, i);

    for(i = 0; i < ((*lowrank)->R).cols(); i++)
        for(j = 0; j < ((*lowrank)->R).rows(); j++)
            r[((*lowrank)->R).rows() * i + j] = ((*lowrank)->R)(j, i);

    return;
}

void initialize_hodlr_c(int N, int M, double eps, HODLR** tree)
{
    (*tree) = new HODLR(N, M, eps);
    return;
}

void assemble_c(HODLR**tree, Kernel** kernel, char* lowrank_method, bool is_sym, bool is_pd)
{
    (*tree)->assemble((*kernel), lowrank_method, is_sym, is_pd);
}

void matmat_product_c(HODLR** tree, double* x, double* b)
{
    Mat b_eig, x_eig;
    x_eig = Eigen::Map<Mat>(x, (*tree)->N, 1);
    b_eig = (*tree)->matmatProduct(x_eig);

    for(int i = 0; i < (*tree)->N; i++)
        b[i] = b_eig(i, 0);

    return;
}

void factorize_c(HODLR** tree)
{
    (*tree)->factorize();
    return;
}

void solve_c(HODLR** tree, double* b, double* x)
{
    Mat b_eig, x_eig;

    b_eig = Eigen::Map<Mat>(b, (*tree)->N, 1);
    x_eig = (*tree)->solve(b_eig);

    for(int i = 0; i < (*tree)->N; i++)
        x[i] = x_eig(i, 0);

    return;
}

void logdeterminant_c(HODLR** tree, double &log_det)
{
    log_det = (*tree)->logDeterminant();
    return;
}

void symm_factor_transpose_prod_c(HODLR** tree, double* x, double* b)
{
    Mat b_eig, x_eig;
    x_eig = Eigen::Map<Mat>(x, (*tree)->N, 1);
    b_eig = (*tree)->symmetricFactorTransposeProduct(x_eig);

    for(int i = 0; i < (*tree)->N; i++)
        b[i] = b_eig(i, 0);

    return;
}

void symm_factor_prod_c(HODLR** tree, double* x, double* b)
{
    Mat b_eig, x_eig;
    x_eig = Eigen::Map<Mat>(x, (*tree)->N, 1);
    b_eig = (*tree)->symmetricFactorProduct(x_eig);

    for(int i = 0; i < (*tree)->N; i++)
        b[i] = b_eig(i, 0);

    return;
}

void get_symm_factor_c(HODLR** tree, double* W_matrix)
{
    Mat temp = (*tree)->getSymmetricFactor();
    // Counter variables:
    int i, j;
    for(i = 0; i < temp.cols(); i++)
        for(j = 0; j < temp.rows(); j++)
            W_matrix[temp.rows() * i + j] = temp(j, i);

    return;
}
