#include "HODLR_Matrix.hpp"
#include "Matrix_Factorizer.hpp"
#include "KDTree.hpp"

// Derived class of HODLR_Matrix which is ultimately
// passed to the HODLR_Tree class:
class Kernel : public HODLR_Matrix 
{

private:
    Vec x, y;

public:

    // Constructor:
    Kernel(int N) : HODLR_Matrix(N) 
    {
        x = 2 * Vec::Ones(N) + Vec::Random(N);
        y = 7 * Vec::Ones(N) + Vec::Random(N);

        std::sort(x.data(), x.data() + x.size(), std::greater<double>());
        std::sort(y.data(), y.data() + y.size());
        std::ofstream myfile;
        myfile.open("x.txt");
        for(int j = 0; j < x.size(); j++)
        {
            myfile << x[j] << std::endl;
        }
        // Closing the file:
        myfile.close();

        myfile.open("y.txt");
        for(int j = 0; j < y.size(); j++)
        {
            myfile << y[j] << std::endl;
        }
        // Closing the file:
        myfile.close();
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        return 1 / abs((x(i)-y(j)));
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    srand(time(NULL));
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
    Kernel* K             = new Kernel(N);
    Matrix_Factorizer* F1 = new Matrix_Factorizer(K, "rookPivoting");
    Matrix_Factorizer* F2 = new Matrix_Factorizer(K, "queenPivoting");
    Matrix_Factorizer* F3 = new Matrix_Factorizer(K, "SVD");
    Matrix_Factorizer* F4 = new Matrix_Factorizer(K, "RRQR");

    Mat B = K->getMatrix(0, 0, N, N);
    Mat L, R, error;

    F1->getFactorization(L, R, tolerance);
    error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization using Rook Pivoting:" << error.cwiseAbs().maxCoeff() << std::endl;

    F2->getFactorization(L, R, tolerance);
    error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization using Queen Pivoting:" << error.cwiseAbs().maxCoeff() << std::endl;

    F3->getFactorization(L, R, tolerance);
    error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization using SVD:" << error.cwiseAbs().maxCoeff() << std::endl;

    F4->getFactorization(L, R, tolerance);
    error = B - L * R.transpose();
    std::cout << "Accuracy of Factorization using RRQR:" << error.cwiseAbs().maxCoeff() << std::endl << std::endl;

    return 0;
}
