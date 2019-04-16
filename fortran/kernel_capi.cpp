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
    K->getMatrixEntry(i, j);
}
