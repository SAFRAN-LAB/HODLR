#ifndef __HODLR_Matrix__
#define __HODLR_Matrix__

#ifdef MKL_ENABLED
    #define EIGEN_USE_MKL_ALL
#endif

#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
// Used to dump data:
#include <fstream>
#include <complex>

using std::cout;
using std::endl;

#ifdef USE_FLOAT
    using dtype=float;
    using dtype_base=float;
    using Mat=Eigen::MatrixXf;
    using Vec=Eigen::VectorXf;
#endif

#ifdef USE_DOUBLE
    using dtype=double;
    using dtype_base=double;
    using Mat=Eigen::MatrixXd;
    using Vec=Eigen::VectorXd;
#endif

#ifdef USE_COMPLEX32
    using dtype=std::complex<float>;
    using dtype_base=float;
    using Mat=Eigen::MatrixXcf;
    using Vec=Eigen::VectorXcf;
    const std::complex<float> I(0.0, 1.0);
#endif

#ifdef USE_COMPLEX64
    using dtype=std::complex<double>;
    using dtype_base=double;
    using Mat=Eigen::MatrixXcd;
    using Vec=Eigen::VectorXcd;
    const std::complex<double> I(0.0, 1.0);
#endif

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
    virtual dtype getMatrixEntry(int j, int k) 
    {
        // FROM EXPERIENCE: Incase the user makes a mistake in 
        // setting the derived class, this warns the user:
        cout << "Returning zero! Ensure that derived class is properly set!" << endl;
        return 0.0;
    }

    Vec getRow(int j, int n_col_start, int n_cols);
    Vec getCol(int k, int n_row_start, int n_rows);
    Mat getMatrix(int j, int k, int n_rows, int n_cols);

    void maxAbsVector(const Vec& v, 
                      const std::set<int>& allowed_indices, 
                      dtype& max, int& index
                     );

    void rookPiv(int n_row_start, int n_col_start, int n_rows, int n_cols, double tolerance, 
                 Mat& L,  Mat& R, int& computed_rank
                );
  
    // Destructor:
    ~HODLR_Matrix() {};
};

#endif /*__HODLR_Matrix__*/
