// This file serves as a gentle introduction to the usage of this library:

#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using std::cout;
using std::endl;

// Derived class of HODLR_Matrix which is ultimately
// passed to the HODLR_Tree class:

class myHODLR_Matrix : public HODLR_Matrix 
{

private:
    VectorXd x;

public:

    // Constructor:
    myHODLR_Matrix(int N) : HODLR_Matrix(N) 
    {
        x = VectorXd::Random(N);
        // This is being sorted to ensure that we get
        // optimal low rank structure:
        std::sort(x.data(),x.data()+x.size());
    };
    
    // In this example, we are illustrating usage using
    // the gaussian kernel:
    double getMatrixEntry(int j, int k) 
    {
        // Value on the diagonal:
        if(j==k) 
        {
            return 10;
        }
        
        // Otherwise:
        else 
        {
            return exp(-(x(j)-x(k))*(x(j)-x(k)));
        }
    }

    // Destructor:
    ~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) 
{
    // Size of the Matrix in consideration:
    int N             = atoi(argv[1]);
    // Declaration of HODLR_Matrix object that abstracts data in Matrix:
    myHODLR_Matrix* A = new myHODLR_Matrix(N);
    // Here it is assumed that size of leaf level is 200
    int n_levels      = log(N/200) / log(2);
    double tolerance  = 1e-12;

    // Variables used in timing:
    double start, end;

    cout << "Fast method..." << endl;
    
    start = omp_get_wtime();
    // Creating a pointer to the HODLR Tree structure:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, A);
    T->assembleTree();
    end   = omp_get_wtime();
    
    cout << "Time for assembly in HODLR form:" << (end - start) << endl;

    // This is used in debugging mainly:
    T->printTreeDetails();

    // Random Matrix to multiply with
    MatrixXd x = MatrixXd::Random(N, 1);
    // Stores the result after multiplication:
    MatrixXd b_fast;
    
    start = omp_get_wtime();
    T->matmatProduct(x, b_fast);
    end   = omp_get_wtime();
    
    cout << "Time for matrix-vector product:" << (end - start) << endl << endl;

    cout << "Exact method..." << endl;

    // What we are doing here is explicitly generating 
    // the matrix from its entries
    start = omp_get_wtime();
    MatrixXd B = A->getMatrix(0,0,N,N);
    end   = omp_get_wtime();
    
    cout << "Time for matrix generation:" << (end-start) << endl;

    start = omp_get_wtime();
    MatrixXd b_exact = B * x;
    end   = omp_get_wtime();
    
    cout << "Time for matrix-vector product:" << (end-start) << endl;
    // Computing the relative error in the solution obtained:
    cout << "Error in the solution is:" << (b_fast-b_exact).norm() / (b_exact.norm()) << endl << endl;

    start = omp_get_wtime();
    T->factorize();
    end   = omp_get_wtime();
    cout << "Time to factorize:" << (end-start) << endl;

    MatrixXd x_fast;
    start  = omp_get_wtime();
    x_fast = T->solve(b_fast);
    end    = omp_get_wtime();
    cout << "Time to solve:" << (end-start) << endl;
    // Computing the relative error:
    cout << "Error in the solution:" << (x_fast - x).norm() / (x.norm()) << endl << endl;

    // Computing log-determinant using Cholesky:
    Eigen::LLT<MatrixXd> P;
    start = omp_get_wtime();
    P.compute(B);
    double log_det = 0.0;
    for(int i=0; i<P.matrixL().rows(); ++i)
    {
        log_det += log(P.matrixL()(i,i));
    }
    end = omp_get_wtime();
    log_det = 2 * log_det;
    cout << "Time to calculate log determinant using Cholesky:" << (end-start) << endl;

    start = omp_get_wtime();
    double log_det_hodlr = T->logDeterminant();
    end = omp_get_wtime();
    cout << "Time to calculate log determinant using HODLR:" << (end-start) << endl;
    cout << "Error in computation:" << fabs(log_det_hodlr - log_det) << endl;
}
