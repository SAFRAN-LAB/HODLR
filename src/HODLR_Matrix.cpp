#include "HODLR_Matrix.hpp"

Eigen::VectorXd HODLR_Matrix::getRow(const int j, const int n_col_start, const int n_cols) 
{
    Eigen::VectorXd row(n_cols);
    #pragma omp parallel for
    for(int k = 0; k < n_cols; k++) 
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

    for(it = allowed_indices.begin(); it != allowed_indices.end(); it++) 
    {
        if(fabs(v(*it))>fabs(max)) 
        {
            index   =   *it;
            max     =   v(index);
        }
    }
}

// This performs ACA upon the input matrix and returns U and V such that
// ||A - UV^T|| < tolerance
void HODLR_Matrix::rookPiv(int n_row_start, int n_col_start, 
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

    double max, gamma, unused_max;

    // Initialize the matrix norm and the the first row index
    double matrix_norm = 0;
    row_ind.push_back(0);
    remaining_row_ind.erase(0);

    // Stores the pivot entry of the considered row / col:
    int pivot;

    // This would get updated:
    computed_rank = 0;
    Eigen::VectorXd row, col;
    // These quantities in finding the stopping criteria:
    double row_squared_norm, row_norm, col_squared_norm, col_norm;

    // So these would be particularly useful for poorly conditioned matrices:
    int max_tries = 20;
    int count;

    // Repeat till the desired tolerance is obtained
    do 
    {
        // Generation of the row
        // Row of the residuum and the pivot column
        // By calling row_ind.back(), we are getting the last pushed number
        row = this->getRow(n_row_start + row_ind.back(), n_col_start, n_cols);
        
        for(int i = 0; i < computed_rank; i++) 
        {
            row = row - u[i](row_ind.back()) * v[i];
        }

        this->maxAbsVector(row, remaining_col_ind, max, pivot);
        count = 0;

        // This randomization is needed if in the middle of the algorithm the 
        // row happens to be exactly the linear combination of the previous rows 
        // upto some tolerance.
        // Alternating upon each call:
        bool eval_at_end = false;
        
        while (fabs(max) < tolerance && 
               count < max_tries   && 
               remaining_col_ind.size() > 0 && 
               remaining_row_ind.size() > 0
              ) 
        {
            row_ind.pop_back();
            int new_row_ind;

            if(eval_at_end)
            {
                new_row_ind = *remaining_row_ind.end();
            }

            else
            {
                new_row_ind = *remaining_row_ind.begin();
            }

            eval_at_end = !eval_at_end;
            row_ind.push_back(new_row_ind);
            remaining_row_ind.erase(new_row_ind);

            // Generation of the row
            // Row of the residuum and the pivot column
            row = this->getRow(n_row_start + new_row_ind, n_col_start, n_cols);
            for(int i = 0; i < computed_rank; i++) 
            {
                row = row - u[i](row_ind.back()) * v[i];
            }

            this->maxAbsVector(row, remaining_col_ind, max, pivot);
            count++;
        }

        // In case it failed to resolve in the previous step, 
        // we break out of the dowhile loop:
        if (count == max_tries || 
            remaining_col_ind.size() == 0 || 
            remaining_row_ind.size() == 0
           )
        {
            break;
        } 

        // Now we will move onto doing this process for the column bases
        // Resetting the count back to zero:
        count = 0;

        col_ind.push_back(pivot);
        remaining_col_ind.erase(pivot);

        // Normalizing constant
        gamma = 1.0 / max;

        // Generation of the column
        // Column of the residuum and the pivot row
        col = this->getCol(n_col_start + col_ind.back(), n_row_start, n_rows);
        
        for(int i = 0; i < computed_rank; i++) 
        {
            col = col - v[i](col_ind.back()) * u[i];
        }

        this->maxAbsVector(col, remaining_row_ind, unused_max, pivot);

        // Repeating the same randomization we carried out for the rows, now for the columns:
        while (fabs(max)<tolerance && 
               count < max_tries && 
               remaining_col_ind.size() >0 && 
               remaining_row_ind.size() >0
              ) 
        {
            col_ind.pop_back();
         
            int new_col_ind;

            if(eval_at_end)
            {
                new_col_ind = *remaining_col_ind.end();
            }

            else
            {
                new_col_ind = *remaining_col_ind.begin();
            }

            eval_at_end = !eval_at_end;

            col_ind.push_back(new_col_ind);
            remaining_col_ind.erase(new_col_ind);

            // Generation of the column
            // Column of the residuum and the pivot row
            col = this->getCol(n_col_start + new_col_ind, n_row_start, n_rows);
            for(int i = 0; i < computed_rank; i++) 
            {
                col = col - v[i](col_ind.back()) * u[i];
            }

            this->maxAbsVector(col, remaining_row_ind, unused_max, pivot);
            count++;
        }

        // In case it failed to resolve in the previous step, 
        // we break out of the dowhile loop:
        if (count == max_tries || 
            remaining_col_ind.size() == 0 || 
            remaining_row_ind.size() == 0
           )
        {
            break;
        } 

        // Resetting back for next iteration:
        count = 0;

        row_ind.push_back(pivot);
        remaining_row_ind.erase(pivot);

        // New vectors
        u.push_back(gamma * col);
        v.push_back(row);

        // New approximation of matrix norm
        row_squared_norm = row.squaredNorm();
        row_norm         = sqrt(row_squared_norm);

        col_squared_norm = col.squaredNorm();
        col_norm         = sqrt(col_squared_norm);

        // Updating the matrix norm:
        matrix_norm = matrix_norm + gamma * gamma * row_squared_norm * col_squared_norm;

        for(int j = 0; j < computed_rank; j++) 
        {
            matrix_norm = matrix_norm + 2.0 * (u[j].dot(u.back())) * (v[j].dot(v.back()));
        }
        
        computed_rank++;
    } 
    while(computed_rank * (n_rows + n_cols) * row_norm * col_norm > 
          fabs(max) * tolerance * matrix_norm && 
          computed_rank < fmin(n_rows, n_cols)
         );

    // If the computed_rank is >= to full-rank
    // then return the trivial full-rank decomposition
    if (computed_rank >= fmin(n_rows, n_cols) - 1) 
    {
        if (n_rows < n_cols) 
        {
            L = Eigen::MatrixXd::Identity(n_rows, n_rows);
            R = getMatrix(n_row_start, n_col_start, n_rows, n_cols).transpose();
            computed_rank = n_rows;
        }

        else 
        {
            L = getMatrix(n_row_start, n_col_start, n_rows, n_cols);
            R = Eigen::MatrixXd::Identity(n_cols, n_cols);
            computed_rank = n_cols;
        }
    }
    
    // This is when ACA has succeeded:
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
