#include "fortran_wrapper.hpp"

void initialize_kernel_object(Kernel** kernel, int N, int dim)
{
    (*kernel) = new Kernel(N, dim);
    return;
}

void get_matrix(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end)
{
    double* temp = ((*kernel)->getMatrix(row_start, col_start, row_end, col_end)).data();
    // HACK: For some reason the first entry of the matrix isn't getting captured:
    matrix[0] = (*kernel)->getMatrixEntry(0, 0);
    for(int i = 1; i < (row_end - row_start) * (col_end - col_start); i++)
        matrix[i] = temp[i];

    return;
}
