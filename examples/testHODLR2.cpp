#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

using std::vector;
using Eigen::MatrixXd;

// Example provided by Michael-Hartmann:
// The matrix is given by:
//       D = Id - M_ij
//  where
//       M_ij = y^(i+j+1) * (i+j)!/(i! j!)
//  with the indices 1 <= i,j <= ldim, and
//       y = 0.5*R/(R+L).
 
double __kernel(double y, int i, int j)
{
    const int l1 = i+1, l2 = j+1;
    return -exp((l1 + l2 + 1) * y + lgamma(1 + l1 + l2) - lgamma(1 + l1) - lgamma(1 + l2));
}

class kernel : public HODLR_Matrix 
{

private:
    double y;

public:

    kernel(unsigned N, double RbyL) : HODLR_Matrix(N) 
    {
        // y = 0.5*R/(R+L) 
        this->y = log(0.5/(1+1./RbyL));
    };

    double getMatrixEntry(int i, int j) 
    {
        // Value on the diagonal:
        if(i == j) 
        {
            return (1 + __kernel(this->y, i, j));
        }
        
        // Otherwise:
        else 
        {
            return __kernel(this->y, i, j);
        }
    };
};

int main(int argc, char *argv[]) 
{
    // exact (analytical) result from
    // Giuseppe Bimonte, Classical Casimir interaction of a perfectly
    // conducting sphere and plate, Phys. Rev. D 95, 065004 (2017)
    // https://doi.org/10.1103/PhysRevD.95.065004
     
    const double exact = -117.00871690447;
    double logdet;

    // Parameters of the problem 
    const double RbyL = 42000; /* R/L */
    const unsigned int dim = 10 * RbyL;

    // Hyper parameters for HODLR 
    const double tolerance = 1e-11;
    const int nLeaf = 100;

    // the "interesting" seed is 1522913970 
    unsigned int seed = 1522913970; // time(NULL);
    if(argc > 1)
        seed = atoi(argv[1]);

    // initialize PRNG 
    printf("initialize PRNG with seed %u\n", seed);
    srand(seed);

    kernel K(dim, RbyL);

    int n_levels  = log(dim / nLeaf) / log(2);
    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, &K);

    // Assemble symmetric matrix 
    T->assembleTree(true);
    // Printing the tree details. Used to debug:
    // T->printTreeDetails();

    // Compute factorization 
    T->factorize();

    // Compute determinant 
    logdet = T->logDeterminant();

    printf("Exact: %.14g\n", exact);
    printf("HODLR: %.14g\n", logdet);
    printf("Relative Error: %g\n", 1 - fabs(logdet / exact));

    return 0;
}
