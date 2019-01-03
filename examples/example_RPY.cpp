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
    int dim;
    double k, T, eta, a;

public:

    // Constructor:
    Kernel(int N, int dim, double k, double T, double eta) : HODLR_Matrix(dim * N) 
    {
        x         = (Mat::Random(N, dim)).real();
        this->dim = dim;
        this->k   = k;
        this->T   = T;
        this->eta = eta;

        // This is being sorted to ensure that we get
        // optimal low rank structure:
        getKDTreeSorted(x, 0);

        // Finding the minimum separation distance:
        double min_r = 1000000;
        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < N; j++)
            {
                double R2 = 0;
                for(int k = 0; k < dim; k++) 
                {
                    R2 += (x(i,k) - x(j,k)) * (x(i,k) - x(j,k));
                }
                
                if(i != j && sqrt(R2) < min_r)
                    min_r = sqrt(R2);
            }
        }

        // We will set a as the minimum of the distances between particles / 2:
        this->a = min_r / 2;
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        if(i / dim == j / dim)
        {
            if(i % dim == j % dim)
                return ((k * T) / (6 * PI * eta * a));
            else
                return 0;
        }

        else 
        {
            dtype R2 = 0;
            Vec r(dim);
            for(int k = 0; k < dim; k++) 
            {
                r(k) = x(i / dim, k) - x(j / dim, k);
                R2  += r(k) * r(k);
            }
            dtype R = sqrt(R2);
            Mat tensor = ((this->k * T) / (8 * PI * eta * R)) * (  Mat::Identity(dim, dim) + r * r.transpose() / R2 
                                                                 + (2 * a * a / (3 * R2)) * (  Mat::Identity(dim, dim) 
                                                                                             - 3 * r * r.transpose() / R2
                                                                                            )
                                                                );
            return tensor(i % dim, j % dim);
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
    // Setting k = T = η = 1
    Kernel* K         = new Kernel(N, dim, 1, 1, 1);
    int n_levels      = log(dim * N / M) / log(2);

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
    Mat x = (Mat::Random(N * dim, 1)).real();
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
    Mat B = K->getMatrix(0, 0, dim * N, dim * N);
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
