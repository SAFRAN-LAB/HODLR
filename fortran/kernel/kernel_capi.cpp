#include "kernel.h"
#include "Kernel.hpp"

KERNEL* create_kernel(int N)
{
    return new KERNEL(N);
}

void delete_kernel(KERNEL* K)
{
    delete K;
}

double get_matrix_entry(const KERNEL* K, int i, int j)
{
    return K->getMatrixEntry(i, j);
}

void get_vector_x(const KERNEL* K, const double* x)
{
    K->getVectorX(x);
}

// double* get_matrix(const KERNEL* K, int row_start, int col_start, int row_end, int col_end)
// {
//     Eigen::MatrixXd a(row_end - row_start, col_end - col_start);

//     A = K->getMatrix(row_start, col_start, row_end, col_end)
//     return a.data()
// }
