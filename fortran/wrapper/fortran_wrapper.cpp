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

void initialize_kernel_object_c(Kernel** kernel, int N)
{
    (*kernel) = new Kernel(N);
    return;
}

void get_matrix_c(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end)
{
    Mat temp = (*kernel)->getMatrix(row_start, col_start, row_end, col_end);
    // Counter variables:
    int i, j;
    for(i = 0; i < temp.cols(); i++)
        for(j = 0; j < temp.rows(); j++)
            matrix[temp.rows() * i + j] = temp(j, i);

    return;
}

void initialize_matrix_factorizer_c(Matrix_Factorizer** factorizer, Kernel** kernel, char* factorization_type)
{   
    removeEndSpaces(factorization_type);
    (*factorizer) = new Matrix_Factorizer((*kernel), factorization_type);
    return;
}

void get_factorization_c(Matrix_Factorizer** factorizer, double* l, double* r, double eps)
{
    Mat U, V;
    (*factorizer)->getFactorization(U, V, eps);

    // Counter variables:
    int i, j;
    for(i = 0; i < U.cols(); i++)
        for(j = 0; j < U.rows(); j++)
            l[U.rows() * i + j] = U(j, i);

    for(i = 0; i < V.rows(); i++)
        for(j = 0; j < V.cols(); j++)
            r[V.cols() * i + j] = V(i, j);

    return;
}

void initialize_hodlr_tree_c(HODLR_Tree** tree, int n_levels, double eps, Matrix_Factorizer** factorizer)
{
    (*tree) = new HODLR_Tree(n_levels, eps, (*factorizer));
    return;
}

void assemble_tree_c(HODLR_Tree** tree, bool is_sym, bool is_pd)
{
    (*tree)->assembleTree(is_sym, is_pd);
    return;
}

void matmat_product_c(HODLR_Tree** tree, double* x, double* b)
{
    Mat b_eig, x_eig;
    x_eig = Eigen::Map<Mat>(x, (*tree)->N, 1);
    b_eig = (*tree)->matmatProduct(x_eig);

    for(int i = 0; i < (*tree)->N; i++)
        b[i] = b_eig(i, 0);

    return;
}

void factorize_c(HODLR_Tree** tree)
{
    (*tree)->factorize();
    return;
}

void solve_c(HODLR_Tree** tree, double* b, double* x)
{
    Mat b_eig, x_eig;

    b_eig = Eigen::Map<Mat>(b, (*tree)->N, 1);
    x_eig = (*tree)->solve(b_eig);

    for(int i = 0; i < (*tree)->N; i++)
        x[i] = x_eig(i, 0);

    return;
}

void logdeterminant_c(HODLR_Tree** tree, double &log_det)
{
    log_det = (*tree)->logDeterminant();
    return;
}
