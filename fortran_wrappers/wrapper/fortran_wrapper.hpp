#ifndef __fortran_wrapper__
#define __fortran_wrapper__

#include "HODLR_Tree.hpp"
#include "Kernel.hpp"

// Extern C needs to be used to prevent name mangling
extern "C"
{
    void initialize_kernel_object(Kernel** kernel, int N, int dim);
    void get_matrix(double* matrix, Kernel** kernel, int row_start, int col_start, int row_end, int col_end);
}

#endif /* (defined(__fortran_wrapper__) */
