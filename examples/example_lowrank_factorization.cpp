#include "HODLR_Matrix.hpp"
#include "Matrix_Factorizer.hpp"
#include "KDTree.hpp"

// Derived class of HODLR_Matrix which is ultimately
// passed to the HODLR_Tree class:
class Kernel : public HODLR_Matrix 
{

private:
    Mat x;

public:

    // Constructor:
    Kernel(int N) : HODLR_Matrix(N) 
    {
        Vec a = Vec::Random(N);
        Vec b = Vec::Random(N);

        x = a * b.transpose();
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        return x(i, j);
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    int N, M, dim;
    double tolerance;

    if(argc < 3)
    {
        std::cout << "All arguments weren't passed to executable!" << std::endl;
        std::cout << "Using Default Arguments:" << std::endl;
        // Size of the Matrix in consideration:
        N          = 1000;
        // Tolerance of problem
        tolerance  = pow(10, -12);
    }

    else
    {
        // Size of the Matrix in consideration:
        N          = atoi(argv[1]);
        // Tolerance of problem
        tolerance  = pow(10, -atoi(argv[2]));
    }

    std::cout << "========================= Problem Parameters =========================" << std::endl;
    std::cout << "Matrix Size :" << N << std::endl;
    std::cout << "Tolerance   :" << tolerance << std::endl << std::endl;

    // Declaration of HODLR_Matrix object that abstracts data in Matrix:
    Kernel* K            = new Kernel(N);
    Matrix_Factorizer* F = new Matrix_Factorizer(K);

    Mat L, R;
    F->getFactorization(L, R, tolerance);

    Mat B     = K->getMatrix(0, 0, N, N);
    Mat error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization:" << error.maxCoeff() << std::endl << std::endl;

    return 0;
}
