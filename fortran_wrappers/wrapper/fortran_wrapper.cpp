#include "fortran_wrapper.hpp"

// This is a helper function:
void removeSpaces(char *str) 
{ 
    // To keep track of non-space character count 
    int count = 0; 
  
    // Traverse the given string. If current character 
    // is not space, then place it at index 'count++' 
    for (int i = 0; str[i]; i++) 
        if (str[i] != ' ') 
            str[count++] = str[i]; // here count is 
                                   // incremented 
    str[count] = '\0'; 
} 

void initialize_kernel_object_c(Kernel** kernel, int N, int dim)
{
    std::cout << "This is the input given to initialize_kernel_object_c" << std::endl;
    std::cout << "kernel = " << kernel << std::endl;
    std::cout << "N = " << N << std::endl;
    std::cout << "dim = " << dim << std::endl;

    (*kernel) = new Kernel(N, dim);
    return;
}

void get_matrix_c(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end)
{

    std::cout << "This is the input given to get_matrix_c" << std::endl;
    std::cout << "kernel = " << kernel << std::endl;
    std::cout << "row_start = " << row_start << std::endl;
    std::cout << "col_start = " << col_start << std::endl;
    std::cout << "row_end = " << row_end << std::endl;
    std::cout << "col_end = " << col_end << std::endl;

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
    std::cout << "This is the input given to initialize_matrix_factorizer_c" << std::endl;
    std::cout << "factorizer = " << factorizer << std::endl;
    std::cout << "kernel = " << kernel << std::endl;
    removeSpaces(factorization_type);
    std::cout << "factorization_type = " << std::string(factorization_type) << " type" << std::endl;
    (*factorizer) = new Matrix_Factorizer((*kernel), std::string(factorization_type));
    return;
}

void get_factorization_c(Matrix_Factorizer** factorizer, double* l, double* r, double eps)
{
    std::cout << "This is the input given to get_factorization_c" << std::endl;
    std::cout << "factorizer = " << factorizer << std::endl;
    std::cout << "l = " << l << std::endl;
    std::cout << "r = " << r << std::endl;
    std::cout << "eps = " << eps << std::endl;

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
    std::cout << "This is the input given to initialize_hodlr_tree_c" << std::endl;
    std::cout << "tree = " << tree << std::endl;
    std::cout << "n_levels = " << n_levels << std::endl;
    std::cout << "eps = " << eps << std::endl;
    std::cout << "factorizer = " << factorizer << std::endl;

    (*tree) = new HODLR_Tree(n_levels, eps, (*factorizer));
    return;
}

void assemble_tree_c(HODLR_Tree** tree, bool is_sym, bool is_pd)
{
    std::cout << "This is the input given to assemble_tree_c" << std::endl;
    std::cout << "tree = " << tree << std::endl;
    std::cout << "is_sym = " << is_sym << std::endl;
    std::cout << "is_pd = " << is_pd << std::endl;

    (*tree)->assembleTree(is_sym, is_pd);
    return;
}

void matmat_product_c(HODLR_Tree** tree, double* x, double* b)
{
    std::cout << "This is the input given to matmat_product_c" << std::endl;
    std::cout << "tree = " << tree << std::endl;

    Mat b_eig, x_eig;
    x_eig = Eigen::Map<Mat>(x, (*tree)->N, 1);
    b_eig = (*tree)->matmatProduct(x_eig);

    for(int i = 0; i < (*tree)->N; i++)
        b[i] = b_eig(i, 0);

    return;
}

// void factorize_c(HODLR_Tree** tree)
// {
//     (*tree)->factorize();
// }

// void solve_c(HODLR_Tree** tree, double* b, double* x)
// {
//     Mat b_eig, x_eig;
//     Eigen::Map<Mat> b_eig()
// }

// void logdeterminant_c(HODLR_Tree** tree, double &log_det)
// {
//     log_det = (*tree)->logDeterminant()
// }
