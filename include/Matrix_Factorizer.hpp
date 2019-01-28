#ifndef __Matrix_Factorizer__
#define __Matrix_Factorizer__

#include <Eigen/Dense>
#include "HODLR_Matrix.hpp"

class Matrix_Factorizer
{
private:
    
    // Used by the rook-pivoting routine:
    void maxAbsVector(const Vec& v, 
                      const std::set<int>& allowed_indices, 
                      dtype& max, int& index
                     );

    void rookPiv(Mat& L,  Mat& R, double rank_or_tolerance,
                 int n_row_start, int n_col_start, 
                 int n_rows, int n_cols
                );

    void SVD(Mat& L,  Mat& R, double rank_or_tolerance,
             int n_row_start, int n_col_start, 
             int n_rows, int n_cols
            );

public:

    HODLR_Matrix* A;
    // Global size of the matrix:
    int N;
    // Type for the factorization:
    std::string type;

    // Constructor:
    Matrix_Factorizer(HODLR_Matrix* A, std::string type = "rookPivoting")
    {
        this->type = type;
        this->A    = A;
        this->N    = A->N;
    }

    void getFactorization(Mat& L,  Mat& R, double rank_or_tolerance = 1e-8,
                          int n_row_start = 0, int n_col_start = 0, 
                          int n_rows = -1, int n_cols = -1
                         );
};

#endif
