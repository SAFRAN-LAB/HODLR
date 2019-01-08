// File used in CI testing

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
            return 100;
        }
        
        // Otherwise:
        else
        {   
            dtype R2 = 0;

            for(int k = 0; k < dim; k++) 
            {
                R2 += (x(i,k) - x(j,k)) * (x(i,k) - x(j,k));
            }

            return exp(-R2);
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

    // Random Matrix to multiply with
    Mat x = (Mat::Random(N, 1)).real();
    // Stores the result after multiplication:
    Mat y_fast, b_fast;

    // Testing fast factorization:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, K);
    bool is_sym = false;
    bool is_pd = false;
    T->assembleTree(is_sym, is_pd);
    T->printTreeDetails();
    T->plotTree();
    
    b_fast      = T->matmatProduct(x);
    Mat B       = K->getMatrix(0, 0, N, N);
    Mat b_exact = B * x;
    assert((b_fast-b_exact).norm() / (b_exact.norm()) < tolerance);

    T->factorize();
    Mat x_fast;
    x_fast = T->solve(b_exact);
    assert((x_fast - x).norm() / (x.norm()) < tolerance);


    dtype log_det;
    Eigen::PartialPivLU<Mat> lu;
    lu.compute(B);
    log_det = 0.0;
    for(int i = 0; i < lu.matrixLU().rows(); i++)
    {
        log_det += log(lu.matrixLU()(i,i));
    }

    dtype log_det_hodlr = T->logDeterminant();
    assert(fabs(1 - fabs(log_det_hodlr/log_det)) < tolerance);
    delete T;

    // Testing fast symmetric factorization:
    T = new HODLR_Tree(n_levels, tolerance, K);
    is_sym = true;
    is_pd = true;
    T->assembleTree(is_sym, is_pd);
    T->printTreeDetails();
    T->plotTree();

    b_fast      = T->matmatProduct(x);
    // Computing the relative error in the solution obtained:
    assert((b_fast-b_exact).norm() / (b_exact.norm()) < tolerance);

    T->factorize();
    x_fast = T->solve(b_exact);
    assert((x_fast - x).norm() / (x.norm()) < tolerance);

    y_fast = T->symmetricFactorTransposeProduct(x);
    b_fast = T->symmetricFactorProduct(y_fast);

    assert((b_fast - b_exact).norm() / (b_exact.norm()) < tolerance);

    Eigen::LLT<Mat> llt;
    llt.compute(B);
    log_det = 0.0;
    for(int i = 0; i < llt.matrixL().rows(); i++)
    {
        log_det += log(llt.matrixL()(i,i));
    }
    log_det *= 2;

    log_det_hodlr = T->logDeterminant();
    assert(fabs(1 - fabs(log_det_hodlr/log_det)) < tolerance);

    // Getting the symmetric factor:
    Mat W  = T->getSymmetricFactor();
    Mat Wt = W.transpose();

    assert((Wt.colPivHouseholderQr().solve(W.colPivHouseholderQr().solve(b_exact)) - x).cwiseAbs().maxCoeff() < tolerance);

    delete T;
    delete K;

    return 0;
}
