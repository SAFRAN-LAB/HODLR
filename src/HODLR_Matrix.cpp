#include "HODLR_Matrix.hpp"

Vec HODLR_Matrix::getRow(const int j, const int n_col_start, const int n_cols) 
{
    Vec row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) 
    {   
        row(k) = this->getMatrixEntry(j,  k + n_col_start);
    }
    
    return row;
}

Vec HODLR_Matrix::getCol(const int k, const int n_row_start, const int n_rows) 
{
    Vec col(n_rows);
    
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) 
    {
        col(j) = this->getMatrixEntry(j + n_row_start, k);
    }
    
    return col;
}


Mat HODLR_Matrix::getMatrix(const int n_row_start, const int n_col_start, 
                            const int n_rows, const int n_cols) 
{
    Mat mat(n_rows, n_cols);
    
    #pragma omp parallel for
    for (int j=0; j < n_rows; ++j) 
    {
        #pragma omp parallel for
        for (int k=0; k < n_cols; ++k) 
        {
            mat(j,k) = this->getMatrixEntry(j + n_row_start, k + n_col_start);
        }
    }

    return mat;
}
