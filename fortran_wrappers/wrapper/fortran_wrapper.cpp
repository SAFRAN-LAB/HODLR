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

void initialize_kernel_object(Kernel** kernel, int N, int dim)
{
    (*kernel) = new Kernel(N, dim);
    std::cout << (*kernel)->getMatrix(0, 0, N, N) << std::endl;
    return;
}

void get_matrix(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end)
{
    std::cout << (*kernel)->getMatrix(row_start, col_start, row_end, col_end) << std::endl;
    double* temp = ((*kernel)->getMatrix(row_start, col_start, row_end, col_end)).data();
    // HACK: For some reason the first entry of the matrix isn't getting captured:
    matrix[0] = (*kernel)->getMatrixEntry(0, 0);
    for(int i = 1; i < (row_end - row_start) * (col_end - col_start); i++)
        matrix[i] = temp[i];

    return;
}

void initialize_matrix_factorizer(Matrix_Factorizer** factorizer, Kernel** kernel, char* factorization_type)
{   
    removeSpaces(factorization_type);
    (*factorizer) = new Matrix_Factorizer((*kernel), std::string(factorization_type));
    return;
}

void get_factorization(Matrix_Factorizer** factorizer, double* l, double* r, double eps)
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

// void initialize_hodlr_tree(HODLR_Tree** tree, int n_levels, double eps, Matrix_Factorizer** factorizer)
// {
//     std::cout << n_levels << std::endl;
//     std::cout << eps << std::endl;


//     (*tree) = new HODLR_Tree(n_levels, eps, (*factorizer));
//     return;
// }
