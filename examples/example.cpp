// Example provided by Michael-Hartmann:
#include "HODLR/HODLR_Matrix.hpp"
#include "HODLR/HODLR.hpp"

#define LMAX 6000
#define GAMMA 0.57721566490153286 // Euler-Mascheroni constant

// Approximation for modified Bessel function I_n(x) when x<<1 (n>=0)
double logIn(int n, double x)
{
    if(n == 0)
        return log1p(x * x /4);
    else
        return (n * log(x / 2) - lgamma(1 + n));
}

// Approximation for modified Bessel function K_n(x) when x<<1 (n>=0)
double logKn(int n, double x)
{
    if(n == 0)
        return (-log(x / 2) - GAMMA);
    else
        return (n * log(2 / x) + lgamma(n) - log(2));
}


class Test_Kernel : public HODLR_Matrix
{
public:
    explicit Test_Kernel (unsigned N) : HODLR_Matrix(N)
    {};

    dtype getMatrixEntry(int i, int j)
    {
        const int m = i - LMAX, n = j - LMAX;

        const double a = 0.01;
        const double b = 0.020020000000000003;

        double term1 = logIn(abs(m), a) - logKn(abs(m), a);
        double term2 = logIn(abs(n), a) - logKn(abs(n), a);
        double term3 = logKn(abs(m + n), b);

        // Value on the diagonal:
        if(i == j)
        {
            return (1-exp(0.5*(term1+term2)+term3));
        }

        // Otherwise:
        else
        {
            return -exp(0.5*(term1+term2)+term3);
        }
    };
};

int main(int argc, char *argv[])
{
    unsigned int seed = time(NULL);
    std::cout << "Seed is: " << seed << std::endl;
    srand(seed);

    // dimension
    int dim = 2 * LMAX + 1;

    // nLeaf is the size of the smallest block at the leaf level
    const unsigned int nLeaf = 100;
    const double tolerance = 1e-12;

    Test_Kernel K(dim);

    std::cout << "========================= Problem Parameters =========================" << std::endl;
    std::cout << "Matrix Size                        :" << dim << std::endl;
    std::cout << "Leaf Size                          :" << nLeaf << std::endl;
    std::cout << "Tolerance                          :" << tolerance << std::endl << std::endl;

    HODLR* T = new HODLR(dim, nLeaf, tolerance);
    T->assemble(&K, "rookPivoting", true, false);

    double logdet_exact = -33.22918044445708; // computed using Python script logdet.py
    double logdet_hodlr = 0;

    // Compute factorization
    T->factorize();
    // Compute log det(Id+M)
    logdet_hodlr = T->logDeterminant();

    std::cout << "Log determinant is: " << std::setprecision(16) << logdet_hodlr << std::endl;
    std::cout << "Log determinant is: " << std::setprecision(16) << logdet_exact << std::endl;
    std::cout << "Relative error is: "  << std::setprecision(16) << fabs(1-logdet_hodlr/logdet_exact) << std::endl;

    return 0;
}
