#include "HODLR_Matrix.hpp"
#include "Eigen/Dense"

Eigen::VectorXd HODLR_Matrix::getRow(const int j, const int n_col_start, const int n_cols) 
{
    Eigen::VectorXd row(n_cols);
    
    #pragma omp parallel for
    for (int k = 0; k < n_cols; k++) 
    {
        row(k) = this->getMatrixEntry(j,  k + n_col_start);
    }
    
    return row;
}

Eigen::VectorXd HODLR_Matrix::getCol(const int k, const int n_row_start, const int n_rows) 
{
    Eigen::VectorXd col(n_rows);
    
    #pragma omp parallel for
    for (int j=0; j<n_rows; ++j) 
    {
        col(j) = this->getMatrixEntry(j + n_row_start, k);
    }
    
    return col;
}


Eigen::MatrixXd HODLR_Matrix::getMatrix(const int n_row_start, const int n_col_start, 
                                        const int n_rows, const int n_cols) 
{
    Eigen::MatrixXd mat(n_rows, n_cols);
    
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


void HODLR_Matrix::maxAbsVector(const Eigen::VectorXd& v, 
                                const std::set<int>& allowed_indices, 
                                double& max, int& index
                               ) 
{
    std::set<int>::iterator it;
    index = *allowed_indices.begin();
    max   = v(index);

    for (it = allowed_indices.begin(); it != allowed_indices.end(); it++) 
    {
        if(fabs(v(*it))>fabs(max)) 
        {
            index   =   *it;
            max     =   v(index);
        }
    }
}


void HODLR_Matrix::rook_Piv(int n_row_start, int n_col_start, 
                            int n_rows, int n_cols, double tolerance, 
                            Eigen::MatrixXd& L, Eigen::MatrixXd& R, int& computed_rank
                           )
{
    // Indices which have been used:
    std::vector<int> row_ind;      
    std::vector<int> col_ind;

    // Indices that are remaining:
    std::set<int> remaining_row_ind;
    std::set<int> remaining_col_ind;
    
    // Bases:
    std::vector<Eigen::VectorXd> u; 
    std::vector<Eigen::VectorXd> v;

    for(int k = 0; k < n_rows; k++) 
    {
        remaining_row_ind.insert(k);
    }
    
    for(int k = 0; k < n_cols; k++) 
    {
        remaining_col_ind.insert(k);
    }

    // Passing the seed as time to the PRNG:
    srand(time(NULL));
    double max, gamma, unused_max;

    // Initialize the matrix norm and the the first row index
    double matrix_norm = 0;
    row_index.push_back(0);
    remaining_row_ind.erase(0);

    int pivot;

    computed_rank = 0;
    Eigen::VectorXd a, row, col;
    double row_squared_norm, row_norm, col_squared_norm, col_norm;

    int max_tries = 10;
    int count;

    // Repeat till the desired tolerance is obtained
    do 
    {
        // Generation of the row
        // Row of the residuum and the pivot column
        // row = A.row(rowIndex.back());
        row = get_Row(n_row_start+rowIndex.back(), n_col_start, n_cols);
        
        for (int l=0; l<computedRank; ++l) 
        {
            row = row - u[l](rowIndex.back()) * v[l];
        }

        pivot   =   max_Abs_Vector(row, remainingColIndex, max);


        count   =   0;

        /// This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows upto some tolerance.
        while (fabs(max)<tolerance && count < max_tries && remainingColIndex.size() >0 && remainingRowIndex.size() >0) {
            rowIndex.pop_back();
            int new_rowIndex    =   *remainingRowIndex.begin();
            rowIndex.push_back(new_rowIndex);
            remainingRowIndex.erase(new_rowIndex);

            /// Generation of the row
            // a    =   A.row(new_rowIndex);
            a   =   get_Row(n_row_start+new_rowIndex, n_col_start, n_cols);
            /// Row of the residuum and the pivot column
            row =   a;
            for (int l=0; l<computedRank; ++l) {
                row =   row-u[l](rowIndex.back())*v[l];
            }
            pivot   =   max_Abs_Vector(row, remainingColIndex, max);
            ++count;
        }

        if (count == max_tries || remainingColIndex.size() == 0 || remainingRowIndex.size() == 0) break;

        count = 0;

        colIndex.push_back(pivot);
        remainingColIndex.erase(pivot);

        /// Normalizing constant
        Gamma   =   1.0/max;

        /// Generation of the column
        // a    =   A.col(colIndex.back());
        a   =   get_Col(n_col_start+colIndex.back(), n_row_start, n_rows);
        /// Column of the residuum and the pivot row
        col =   a;
        for (int l=0; l<computedRank; ++l) {
            col =   col-v[l](colIndex.back())*u[l];
        }
        pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);

        /// This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
        while (fabs(max)<tolerance && count < max_tries && remainingColIndex.size() >0 && remainingRowIndex.size() >0) {
            colIndex.pop_back();
            int new_colIndex    =   *remainingColIndex.begin();
            colIndex.push_back(new_colIndex);
            remainingColIndex.erase(new_colIndex);

            /// Generation of the column
            // a    =   A.col(new_colIndex);
            a   =   get_Col(n_col_start+new_colIndex, n_row_start, n_rows);

            /// Column of the residuum and the pivot row
            col =   a;
            for (int l=0; l<computedRank; ++l) {
                col =   col-u[l](colIndex.back())*v[l];
            }
            pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);
            ++count;
            std::cout << count << "\n";
        }

        if (count == max_tries || remainingColIndex.size() == 0 || remainingRowIndex.size() == 0) break;

        count = 0;

        rowIndex.push_back(pivot);
        remainingRowIndex.erase(pivot);

        /// New vectors
        u.push_back(Gamma*col);
        v.push_back(row);

        /// New approximation of matrix norm
        row_Squared_Norm    =   row.squaredNorm();
        row_Norm            =   sqrt(row_Squared_Norm);

        col_Squared_Norm    =   col.squaredNorm();
        col_Norm            =   sqrt(col_Squared_Norm);

        matrix_Norm         =   matrix_Norm +   Gamma*Gamma*row_Squared_Norm*col_Squared_Norm;

        for (int j=0; j<computedRank; ++j) {
            matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
        }
        ++computedRank;
    } 
    while(computed_rank * (n_rows + n_cols) * row_norm * col_norm > 
          fabs(max) * tolerance * matrix_norm && 
          computed_rank < fmin(n_rows, n_cols)
         );

    // If the computedRank is >= to full-rank
    // then return the trivial full-rank decomposition
    if (computed_rank >= fmin(n_rows, n_cols) - 1) 
    {
        if (n_rows < n_cols) 
        {
            L = Eigen::MatrixXd::Identity(n_rows, n_rows);
            R = get_Matrix(n_row_start, n_col_start, n_rows, n_cols).transpose();
            computed_rank = n_rows;
        }

        else 
        {
            L = get_Matrix(n_row_start, n_col_start, n_rows, n_cols);
            R = Eigen::MatrixXd::Identity(n_cols, n_cols);
            computed_rank = n_cols;
        }
    }
    
    else 
    {
        L = Eigen::MatrixXd(n_rows, computed_rank);
        R = Eigen::MatrixXd(n_cols, computed_rank);
        
        for (int j = 0; j < computed_rank; j++) 
        {
            L.col(j) = u[j];
            R.col(j) = v[j];
        }
    }
}
