#ifndef __HODLR_Matrix__
#define __HODLR_Matrix__

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>

using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::cout;
using std::endl;

// Default argument for certain functions:
static VectorXd EMPTY_VECTOR = VectorXd();

class HODLR_Matrix 
{

// By declaring HODLR_Tree as a friend, the object
// can access private members of this particular class
friend class HODLR_Tree;

private:
    // Size of the matrix:
    int N;

public:

    // Constructor:
    HODLR_Matrix(int N)
    {
        this->N = N;
    }

    // Returns individual matrix 
    virtual double getMatrixEntry(int j, int k) 
    {
        cout << "Returning zero! Ensure that derived class is properly set!" << endl;
        return 0.0;
    }

    Eigen::VectorXd getRow(int j, int n_col_start, int n_cols);
    Eigen::VectorXd getCol(int k, int n_row_start, int n_rows);
    Eigen::MatrixXd getMatrix(int j, int k, int n_rows, int n_cols);

    void maxAbsVector(const Eigen::VectorXd& v, 
                      const std::set<int>& allowed_indices, 
                      double& max, int& index
                     );
  
    void rookPiv(int n_row_start, int n_col_start, int n_rows, int n_cols, double tolerance, 
                 Eigen::MatrixXd& L, Eigen::MatrixXd& R, int& computed_rank
                );
  
    // Destructor:
    ~HODLR_Matrix() {};
};

#endif /*__HODLR_Matrix__*/
