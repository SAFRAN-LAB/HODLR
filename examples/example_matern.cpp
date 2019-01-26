#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

// Taking Matern Kernel for p = 2:
// K(r) = σ^2 * (1 + sqrt(5) * r / ρ + 5/3 * (r / ρ)^2) * exp(-sqrt(5) * r / ρ)
class Kernel : public HODLR_Matrix 
{

private:
    Mat x;
    double sigma_squared, rho;

public:

    // Constructor:
    Kernel(int N, double sigma, double rho) : HODLR_Matrix(N) 
    {
        x                   = (Vec::Random(N)).real();
        this->sigma_squared = sigma * sigma;
        this->rho           = rho;
        // This is being sorted to ensure that we get
        // optimal low rank structure:
        std::sort(x.data(), x.data() + x.size());
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        double R_by_rho = fabs(x(i) - x(j)) / rho;        
        return sigma_squared * (1 + sqrt(5) * R_by_rho + 5/3 * (R_by_rho * R_by_rho)) * exp(-sqrt(5) * R_by_rho);
    }

    // Destructor:
    ~Kernel() {};
};

int main(int argc, char* argv[]) 
{
    int N, M;
    double tolerance;

    if(argc < 4)
    {
        std::cout << "All arguments weren't passed to executable!" << std::endl;
        std::cout << "Using Default Arguments:" << std::endl;
        // Size of the Matrix in consideration:
        N          = 10000;
        // Size of Matrices at leaf level:
        M          = 200;
        // Tolerance of problem
        tolerance  = pow(10, -12);
    }

    else
    {
        // Size of the Matrix in consideration:
        N          = atoi(argv[1]);
        // Size of Matrices at leaf level:
        M          = atoi(argv[2]);
        // Tolerance of problem
        tolerance  = pow(10, -atoi(argv[3]));
    }

    // Declaration of HODLR_Matrix object that abstracts data in Matrix:
    // Taking σ = 10, ρ = 5:
    Kernel* K            = new Kernel(N, 10, 5);
    Matrix_Factorizer* F = new Matrix_Factorizer(K);
    int n_levels         = log(N / M) / log(2);

    std::cout << "========================= Problem Parameters =========================" << std::endl;
    std::cout << "Matrix Size                        :" << N << std::endl;
    std::cout << "Leaf Size                          :" << M << std::endl;
    std::cout << "Number of Levels in Tree           :" << n_levels << std::endl;
    std::cout << "Tolerance                          :" << tolerance << std::endl << std::endl;

    // Variables used in timing:
    double start, end;

    std::cout << "Fast method..." << std::endl;
    
    start = omp_get_wtime();
    // Creating a pointer to the HODLR Tree structure:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, F);
    bool is_sym = true;
    bool is_pd  = true;
    T->assembleTree(is_sym, is_pd);
    end = omp_get_wtime();
    std::cout << "Time for assembly in HODLR form:" << (end - start) << std::endl;

    // Random Matrix to multiply with
    Mat x = (Mat::Random(N, 1)).real();
    // Stores the result after multiplication:
    Mat y_fast, b_fast;
    
    start  = omp_get_wtime();
    b_fast = T->matmatProduct(x);
    end    = omp_get_wtime();
    
    std::cout << "Time for matrix-vector product:" << (end - start) << std::endl << std::endl;

    start = omp_get_wtime();
    T->factorize();
    end   = omp_get_wtime();
    std::cout << "Time to factorize:" << (end-start) << std::endl;

    Mat x_fast;
    start  = omp_get_wtime();
    x_fast = T->solve(b_fast);
    end    = omp_get_wtime();

    std::cout << "Time to solve:" << (end-start) << std::endl;

    if(is_sym == true && is_pd == true)
    {
        start  = omp_get_wtime();
        y_fast = T->symmetricFactorTransposeProduct(x);
        end    = omp_get_wtime();
        std::cout << "Time to calculate product of factor transpose with given vector:" << (end - start) << std::endl;
        
        start  = omp_get_wtime();
        b_fast = T->symmetricFactorProduct(y_fast);
        end    = omp_get_wtime();
        std::cout << "Time to calculate product of factor with given vector:" << (end - start) << std::endl;        
    }
        
    start = omp_get_wtime();
    dtype log_det_hodlr = T->logDeterminant();
    end = omp_get_wtime();
    std::cout << "Time to calculate log determinant using HODLR:" << (end-start) << std::endl;

    // Direct method:
    start = omp_get_wtime();
    Mat B = K->getMatrix(0, 0, N, N);
    end   = omp_get_wtime();

    if(is_sym == true && is_pd == true)
    {
        start = omp_get_wtime();
        Eigen::LLT<Mat> llt;
        llt.compute(B);
        end = omp_get_wtime();
        std::cout << "Time to calculate LLT Factorization:" << (end-start) << std::endl;
    }

    else
    {
        start = omp_get_wtime();
        Eigen::PartialPivLU<Mat> lu;
        lu.compute(B);
        end = omp_get_wtime();
        std::cout << "Time to calculate LU Factorization:" << (end-start) << std::endl;        
    }

    delete K;
    delete T;

    return 0;
}
