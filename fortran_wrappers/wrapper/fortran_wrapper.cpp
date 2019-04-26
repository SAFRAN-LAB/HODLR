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
    return;
}

void get_matrix(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end)
{
    Mat temp = (*kernel)->getMatrix(row_start, col_start, row_end, col_end);
    // Counter variables:
    int i, j;
    for(i = 0; i < temp.cols(); i++)
        for(j = 0; j < temp.rows(); j++)
            matrix[temp.rows() * i + j] = temp(j, i);

    return;
}

void initialize_matrix_factorizer(Matrix_Factorizer** factorizer, Kernel** kernel, char* factorization_type)
{   
    removeSpaces(factorization_type);
    std::cout << "The factorization type selected is " << factorization_type << " type" << std::endl;
    (*factorizer) = new Matrix_Factorizer((*kernel), std::string(factorization_type));
    return;
}

void get_factorization(Matrix_Factorizer** factorizer, double* l, double* r, int eps)
{
    std::cout << eps << std::endl;
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

void initialize_hodlr_tree(HODLR_Tree** tree, int n_levels, double eps, Matrix_Factorizer** factorizer)
{
    std::cout << n_levels << std::endl;
    std::cout << eps << std::endl;

    (*tree) = new HODLR_Tree(n_levels, eps, (*factorizer));
    return;
}
