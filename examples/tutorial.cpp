// This file serves as a gentle introduction to the usage of this library:
// It is also used in CI testing:

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

// Derived class of HODLR_Matrix which is ultimately
// passed to the HODLR_Tree class:
class Kernel : public HODLR_Matrix 
{

private:
    Mat x;

public:

    // Constructor:
    Kernel(int N, int dim) : HODLR_Matrix(N) 
    {
        x = (Mat::Random(N, dim)).real();
        // This is being sorted to ensure that we get
        // optimal low rank structure:
        getKDTreeSorted(x, 0);
    };
    
    // In this example, we are illustrating usage using
    // the gaussian kernel:
    dtype getMatrixEntry(int i, int j) 
    {
        size_t dim = x.cols();

        // Value on the diagonal:
        if(i == j)
        {
            return 10;
        }
        
        // Otherwise:
        else
        {   
            dtype R, R2;
            // Initializing:
            R = R2 = 0;

            for(int k = 0; k < dim; k++) 
            {
                R2 += (x(i,k) - x(j,k)) * (x(i,k) - x(j,k));
            }

            R = sqrt(R2);

            // Exponential: exp(-R)
            // return exp(-R);
            // Gaussian Kernel: e(-R^2)
            // return exp(-R2);
            // Sinc Kernel: sin(R) / R
            // return (sin(R) / R);
            // // Quadric Kernel: (1 + R^2)
            // return (1 + R2);
            // // Inverse-Quadric Kernel: 1 / (1 + R^2)
            // return (1 / (1 + R2));
            // // Multi-Quadric Kernel: sqrt(1 + R^2)
            // return sqrt(1 + R2);
            // // Inverse-Multiquadric Kernel: 1 / sqrt(1 + R^2)
            // return (1 / sqrt(1 + R2));
            // // Log(R) Kernel:
            // return log(R);
            // // R^2 log(R) Kernel:
            return (R2 * log(R));
            // // 1 / R Kernel:
            // return (1 / R);
            // // log(1 + R) Kernel:
            // return log(1 + R);
        }
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    // Size of the Matrix in consideration:
    int N             = atoi(argv[1]);
    // Size of Matrices at leaf level:
    int M             = atoi(argv[2]);
    // Dimensionality of the problem:
    int dim           = atoi(argv[3]);
    // Tolerance of problem
    double tolerance  = pow(10, -atoi(argv[4]));
    // Declaration of HODLR_Matrix object that abstracts data in Matrix:
    Kernel* K         = new Kernel(N, dim);
    int n_levels      = log(N / M) / log(2);

    // Variables used in timing:
    double start, end;

    cout << "Fast method..." << endl;
    
    start = omp_get_wtime();
    // Creating a pointer to the HODLR Tree structure:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, K);
    // If we are assembling a symmetric matrix:
    bool is_sym = true;
    // If we know that the matrix is also PD:
    // By toggling this flag to true, the factorizations are performed using Cholesky
    // Useful when you want the factorization as WW^T 
    bool is_pd = true;
    T->assembleTree(is_sym, is_pd);
    end = omp_get_wtime();
    cout << "Time for assembly in HODLR form:" << (end - start) << endl;

    // These are used in debugging mainly:
    // T->plotTree();
    // T->printTreeDetails();

    // Random Matrix to multiply with
    Mat x = (Mat::Random(N, 1)).real();
    // Stores the result after multiplication:
    Mat y_fast, b_fast;
    
    start  = omp_get_wtime();
    b_fast = T->matmatProduct(x);
    end    = omp_get_wtime();
    
    cout << "Time for matrix-vector product:" << (end - start) << endl << endl;
    // cout << "Exact method..." << endl;

    // What we are doing here is explicitly generating 
    // the matrix from its entries
    // start = omp_get_wtime();
    // Mat B = K->getMatrix(0, 0, N, N);
    // end   = omp_get_wtime();
    
    // cout << "Time for matrix generation:" << (end-start) << endl;

    // start = omp_get_wtime();
    // Mat b_exact = B * x;
    // end   = omp_get_wtime();
    
    // cout << "Time for matrix-vector product:" << (end-start) << endl;
    // // Computing the relative error in the solution obtained:
    // cout << "Error in the solution is:" << (b_fast-b_exact).norm() / (b_exact.norm()) << endl << endl;

    // assert((b_fast-b_exact).norm() / (b_exact.norm()) < tolerance);

    start = omp_get_wtime();
    T->factorize();
    end   = omp_get_wtime();
    cout << "Time to factorize:" << (end-start) << endl;

    Mat x_fast;
    start  = omp_get_wtime();
    x_fast = T->solve(b_fast);
    end    = omp_get_wtime();

    cout << "Time to solve:" << (end-start) << endl;
    // Computing the relative error:
    // cout << "Error in the solution:" << (x_fast - x).norm() / (x.norm()) << endl << endl;

    // assert((x_fast - x).norm() / (x.norm()) < tolerance);

    // Checking symmetric factor product:
    if(is_sym == true && is_pd == true)
    {
        start  = omp_get_wtime();
        y_fast = T->symmetricFactorTransposeProduct(x);
        end    = omp_get_wtime();
        cout << "Time to calculate product of factor transpose with given vector:" << (end - start) << endl;
        
        start  = omp_get_wtime();
        b_fast = T->symmetricFactorProduct(y_fast);
        end    = omp_get_wtime();
        cout << "Time to calculate product of factor with given vector:" << (end - start) << endl;
        
        // cout << "Error in the solution is:" << (b_fast - b_exact).norm() / (b_exact.norm()) << endl << endl;
    }

    // assert((b_fast - b_exact).norm() / (b_exact.norm()) < tolerance);

    // dtype log_det;
    // // Computing log-determinant using Cholesky:
    // if(is_sym == true && is_pd == true)
    // {
    //     Eigen::LLT<Mat> llt;
    //     start = omp_get_wtime();
    //     llt.compute(B);
    //     log_det = 0.0;
    //     for(int i = 0; i < llt.matrixL().rows(); i++)
    //     {
    //         log_det += log(llt.matrixL()(i,i));
    //     }
    //     log_det *= 2;
    //     end = omp_get_wtime();
    //     cout << "Time to calculate log determinant using Cholesky:" << (end - start) << endl;
    //     cout << "Calculated Log Determinant:" << log_det << endl;
    // }

    // // Computing log-determinant using LU:
    // else
    // {
    //     Eigen::PartialPivLU<Mat> lu;
    //     start = omp_get_wtime();
    //     lu.compute(B);
    //     log_det = 0.0;
    //     for(int i = 0; i < lu.matrixLU().rows(); i++)
    //     {
    //         log_det += log(lu.matrixLU()(i,i));
    //     }
    //     end = omp_get_wtime();
    //     cout << "Time to calculate log determinant using LU:" << (end - start) << endl;
    //     cout << "Calculated Log Determinant:" << log_det << endl;
    // }

    start = omp_get_wtime();
    dtype log_det_hodlr = T->logDeterminant();
    end = omp_get_wtime();
    cout << "Time to calculate log determinant using HODLR:" << (end-start) << endl;
    // cout << "Calculated Log Determinant:" << log_det_hodlr << endl;
    // cout << "Relative Error in computation:" << fabs(1 - fabs(log_det_hodlr/log_det)) << endl;

    // assert(fabs(1 - fabs(log_det_hodlr/log_det)) < tolerance);

    // // Getting the symmetric factor:
    // if(is_sym == true && is_pd == true)
    // {
    //     Mat W  = T->getSymmetricFactor();
    //     Mat Wt = W.transpose();

    //     assert((Wt.colPivHouseholderQr().solve(W.colPivHouseholderQr().solve(b_exact)) - x).cwiseAbs().maxCoeff() < tolerance);
    // }

    delete K;
    delete T;

    return 0;
}
