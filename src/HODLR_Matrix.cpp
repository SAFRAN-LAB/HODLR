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

Vec HODLR_Matrix::getDiag1(const int n_row_start, const int n_col_start, 
                           const int n_rows, const int n_cols) 
{
    int N = std::max(n_rows, n_cols);
    Vec diag(N);

    int row_ind, col_ind;
    #pragma omp parallel for
    for (int j = 0; j < N; ++j) 
    {   
        if(n_cols > n_rows)
        {
            row_ind = this->mod(n_row_start - n_col_start + j, n_rows);
            col_ind = j;
        }

        else
        {
            row_ind = j;
            col_ind = this->mod(n_col_start - n_row_start + j, n_cols);
        }
        
        diag(j) = this->getMatrixEntry(row_ind, col_ind);
    }
    
    return diag;
}

Vec HODLR_Matrix::getDiag2(const int n_row_start, const int n_col_start, 
                           const int n_rows, const int n_cols) 
{
    int N = std::max(n_rows, n_cols);
    Vec diag(N);
    
    int row_ind, col_ind;
    #pragma omp parallel for
    for (int j = 0; j < N; ++j) 
    {   
        if(n_cols > n_rows)
        {
            row_ind = this->mod(n_row_start + n_col_start - j, n_rows);
            col_ind = j;
        }

        else
        {
            row_ind = j;
            col_ind = this->mod(n_col_start + n_row_start - j, n_cols);
        }

        diag(j) = this->getMatrixEntry(row_ind, col_ind);
    }
    
    return diag;
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
