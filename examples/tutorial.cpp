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
            return exp(-R2);
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
            // return (R2 * log(R));
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

    // Storing Time Taken:
    double hodlr_time, exact_time;
    cout << "========================= Assembly Time =========================" << endl;
    start = omp_get_wtime();
    // Creating a pointer to the HODLR Tree structure:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, K);
    // If we are assembling a symmetric matrix:
    bool is_sym = true;
    // If we know that the matrix is also PD:
    // By setting the matrix to be symmetric-positive definite, 
    // we trigger the fast symmetric factorization method to be used
    // In all other cases the fast factorization method is used
    bool is_pd = true;
    T->assembleTree(is_sym, is_pd);
    end = omp_get_wtime();
    hodlr_time = (end - start);
    cout << "Time for assembly in HODLR form    :" << hodlr_time << endl;

    // What we are doing here is explicitly generating 
    // the matrix from its entries
    start = omp_get_wtime();
    Mat B = K->getMatrix(0, 0, N, N);
    end   = omp_get_wtime();
    exact_time = (end - start);
    cout << "Time for direct matrix generation  :" << exact_time << endl;
    cout << "Magnitude of Speed-Up              :" << (exact_time / hodlr_time) << endl << endl;

    // These are mainly used in development and debugging:
    // Used to visualize the rank structure of the considered kernel:
    // T->plotTree();
    // Prints the details of all the nodes in the tree:
    // T->printTreeDetails();

    // Random Matrix to multiply with
    Mat x = (Mat::Random(N, 1)).real();
    // Stores the result after multiplication:
    Mat y_fast, b_fast;

    cout << "========================= Matrix-Vector Multiplication =========================" << endl;
    start  = omp_get_wtime();
    b_fast = T->matmatProduct(x);
    end    = omp_get_wtime();
    hodlr_time = (end - start);
    cout << "Time for MatVec in HODLR form      :" << hodlr_time << endl;

    start = omp_get_wtime();
    Mat b_exact = B * x;
    end   = omp_get_wtime();
    exact_time = (end - start);
    cout << "Time for direct MatVec             :" << exact_time << endl;
    cout << "Magnitude of Speed-Up              :" << (exact_time / hodlr_time) << endl;
    // Computing the relative error in the solution obtained:
    cout << "Error in the solution is           :" << (b_fast-b_exact).norm() / (b_exact.norm()) << endl << endl;

    cout << "========================= Factorization =========================" << endl;
    start = omp_get_wtime();
    T->factorize();
    end   = omp_get_wtime();
    hodlr_time = (end - start);
    cout << "Time to factorize HODLR form       :" << hodlr_time << endl;

    Eigen::LLT<Mat> llt;
    Eigen::PartialPivLU<Mat> lu;
    
    // Factorizing using Cholesky:
    if(is_sym == true && is_pd == true)
    {
        start = omp_get_wtime();
        llt.compute(B);
        end = omp_get_wtime();
        exact_time = (end - start);
        cout << "Time for Cholesky Factorization    :" << exact_time << endl;
    }

    // Factorizing using LU:
    else
    {
        start = omp_get_wtime();
        lu.compute(B);
        end = omp_get_wtime();
        exact_time = (end - start);
        cout << "Time for LU Factorization          :" << exact_time << endl;
    }

    cout << "Magnitude of Speed-Up              :" << (exact_time / hodlr_time) << endl << endl;

    cout << "========================= Solving =========================" << endl;
    Mat x_fast;
    start  = omp_get_wtime();
    x_fast = T->solve(b_exact);
    end    = omp_get_wtime();
    hodlr_time = (end - start);
    cout << "Time to solve HODLR form           :" << hodlr_time << endl;

    // This is intentionally declared to be different from x:
    Mat x_exact;
    if(is_sym == true && is_pd == true)
    {
        start   = omp_get_wtime();
        x_exact = llt.solve(b_exact);
        end     = omp_get_wtime();
    }

    else
    {
        start   = omp_get_wtime();
        x_exact = lu.solve(b_exact);
        end     = omp_get_wtime();
    }

    exact_time = (end - start);
    cout << "Time taken to solve exactly        :" << exact_time << endl;
    cout << "Magnitude of Speed-Up              :" << (exact_time / hodlr_time) << endl;
    // Computing the relative error:
    cout << "Forward Error(HODLR)               :" << (x_fast - x).norm() / (x.norm()) << endl;
    cout << "Forward Error(Exact)               :" << (x_exact - x).norm() / (x.norm()) << endl;
    cout << "Backward Error(HODLR)              :" << (T->matmatProduct(x_fast) - b_exact).norm() / b_exact.norm() << endl;
    cout << "Backward Error(Exact)              :" << (B * x_exact - b_exact).norm() / b_exact.norm() << endl << endl;

    cout << "========================= Determinant Computation =========================" << endl;
    // Gets the log(determinant) of the matrix abstracted through Kernel:
    start = omp_get_wtime();
    dtype log_det_hodlr = T->logDeterminant();
    end = omp_get_wtime();
    cout << "Time taken for HODLR form          :" << (end-start) << endl;
    cout << "Calculated Log Determinant(HODLR)  :" << log_det_hodlr << endl;

    dtype log_det;
    // Computing log-determinant using Cholesky:
    if(is_sym == true && is_pd == true)
    {
        log_det = 0.0;
        for(int i = 0; i < llt.matrixL().rows(); i++)
        {
            log_det += log(fabs(llt.matrixL()(i,i)));
        }
        log_det *= 2;
        cout << "Calculated Log Determinant(Exact)  :" << log_det << endl;
    }

    // Computing log-determinant using LU:
    else
    {
        log_det = 0.0;
        for(int i = 0; i < lu.matrixLU().rows(); i++)
        {
            log_det += log(fabs(lu.matrixLU()(i,i)));
        }
        cout << "Calculated Log Determinant(Exact)  :" << log_det << endl;
    }

    cout << "Relative Error in computation      :" << fabs(1 - fabs(log_det_hodlr/log_det)) << endl << endl;

    // These routines are specific to the fast symmetric factorization routine:
    // Checking symmetric factor product:
    if(is_sym == true && is_pd == true)
    {
        cout << "========================= Multiplication With Symmetric Factor =========================" << endl;
        // We set y = W^T x
        start  = omp_get_wtime();
        y_fast = T->symmetricFactorTransposeProduct(x);
        end    = omp_get_wtime();
        cout << "Time to calculate product of factor transpose with given vector:" << (end - start) << endl;
        
        // b = W y = W W^T x = B * x
        start  = omp_get_wtime();
        b_fast = T->symmetricFactorProduct(y_fast);
        end    = omp_get_wtime();
        cout << "Time to calculate product of factor with given vector          :" << (end - start) << endl;

        cout << "Error in the solution is                                       :" << (b_fast - b_exact).norm() / (b_exact.norm()) << endl << endl;
    }

    // If we want to explicitly build the symmetric factor matrix, then we can call this command
    if(is_sym == true && is_pd == true)
    {
        Mat W = T->getSymmetricFactor();
    }

    delete K;
    delete T;

    return 0;
}
