#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

#define PI 3.141592653589793238462643383
// RPY Tensor:
// K(r) = ((kT) / (6πηa))     * ((1 - 9 * ||r|| / (32a))I + 3 / (32a) (r ⊗ r) / ||r|| if ||r|| <  2a
// K(r) = ((kT) / (8πη||r||)) * (I + (r ⊗ r) / ||r|| + (2a**2 / 3||r||**2) * (I - 3 (r ⊗ r) / ||r||))  if ||r|| >= 2a

class Kernel : public HODLR_Matrix 
{

private:
    Mat x;
    double k, T, eta, a;

public:

    // Constructor:
    Kernel(int N, int dim, double k, double T, double eta, double a) : HODLR_Matrix(dim * N) 
    {
        x         = (Mat::Random(N, dim)).real();
        this->k   = k;
        this->T   = T;
        this->eta = eta;
        this->a   = a;

        // This is being sorted to ensure that we get
        // optimal low rank structure:
        getKDTreeSorted(x, 0);
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        dtype normR, normR2;
        // Initializing:
        normR = normR2 = 0;

        for(int k = 0; k < dim; k++) 
        {
            normR2 += (x(i,k) - x(j,k)) * (x(i,k) - x(j,k));
        }

        normR = sqrt(normR2);

        if(normR > 2 * a)
        {
            if(i % dim == 0)
            {
                if(j % dim == 0)
                    return 
            }
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
    // Tolerance of problem
    int dim           = atoi(argv[3]);
    double tolerance  = pow(10, -atoi(argv[4]));
    // Declaration of HODLR_Matrix object that abstracts data in Matrix:
    // Taking σ = 10, ρ = 5:
    Kernel* K         = new Kernel(N, dim, 10, 5);
    int n_levels      = log(N / M) / log(2);

    // Variables used in timing:
    double start, end;

    cout << "Fast method..." << endl;
    
    start = omp_get_wtime();
    // Creating a pointer to the HODLR Tree structure:
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, K);
    bool is_sym = true;
    bool is_pd  = false;
    T->assembleTree(is_sym, is_pd);
    end = omp_get_wtime();
    cout << "Time for assembly in HODLR form:" << (end - start) << endl;

    // Random Matrix to multiply with
    Mat x = (Mat::Random(N, 1)).real();
    // Stores the result after multiplication:
    Mat y_fast, b_fast;
    
    start  = omp_get_wtime();
    b_fast = T->matmatProduct(x);
    end    = omp_get_wtime();
    
    cout << "Time for matrix-vector product:" << (end - start) << endl << endl;

    start = omp_get_wtime();
    T->factorize();
    end   = omp_get_wtime();
    cout << "Time to factorize:" << (end-start) << endl;

    Mat x_fast;
    start  = omp_get_wtime();
    x_fast = T->solve(b_fast);
    end    = omp_get_wtime();

    cout << "Time to solve:" << (end-start) << endl;

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
    }
        
    start = omp_get_wtime();
    dtype log_det_hodlr = T->logDeterminant();
    end = omp_get_wtime();
    cout << "Time to calculate log determinant using HODLR:" << (end-start) << endl;

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
        cout << "Time to calculate LLT Factorization:" << (end-start) << endl;
    }

    else
    {
        start = omp_get_wtime();
        Eigen::PartialPivLU<Mat> lu;
        lu.compute(B);
        end = omp_get_wtime();
        cout << "Time to calculate LU Factorization:" << (end-start) << endl;        
    }

    delete K;
    delete T;

    return 0;
}
