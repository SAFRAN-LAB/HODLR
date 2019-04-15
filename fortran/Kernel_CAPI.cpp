#include "kernel_wrap.h"
#include "Kernel.hpp"

#include <iostream>

KERNEL* create_kernel(int N)
{
    return new Kernel(N);
}

void delete_kernel(KERNEL* K)
{
    delete K;
}

double* get_matrix(const KERNEL* K, int i, int j, int n_rows, int n_cols)
{
    Mat A = K->getMatrix(i, j, n_rows, n_cols);
    return A.data();
}
