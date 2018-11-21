#include <iostream>
#include <stdio.h>
#include <vector>
#include <Eigen/Dense>
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
 
double __kernel(double y, unsigned int i, unsigned int j)
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

    double getMatrixEntry(const unsigned i, const unsigned j) 
    {
        return __kernel(this->y, i, j);
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
    const double RbyL = 4200; /* R/L */
    const unsigned int dim = 10 * RbyL;

    // Hyper parameters for HODLR 
    const double tolerance = 1e-11;
    const int nLeaf = 100;

    // the "interesting" seed is 1522913970 
    unsigned int seed = 1522913970;// time(NULL);
    if(argc > 1)
        seed = atoi(argv[1]);

    /* initialize PRNG */
    printf("initialize PRNG with seed %u\n", seed);
    srand(seed);

    kernel K(dim, RbyL);

    // int n_levels  = log(dim / nLeaf) / log(2);
    // HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, &K);

    // // Random Matrix to multiply with
    MatrixXd x = MatrixXd::Random(dim, 1);
    // // Stores the result after multiplication:
    // MatrixXd b_fast;

    // T->matmatProduct(x, b_fast);
    MatrixXd B = K.getMatrix(0,0,dim,dim);
    MatrixXd b_exact = B * x;
    // std::cout << "Error in the solution is:" << (b_fast-b_exact).norm() / (b_exact.norm()) << std::endl << std::endl;

    // T->assembleTree();
    // T->factorize();
    // logdet = T->logDeterminant();
    
    /* compute diagonal entries */
    // const double y = log(0.5/(1+1./RbyL));
    // VectorXd diagonal = VectorXd::Ones(dim);
    // for(int n = 0; n < dim; n++)
    //     diagonal(n) = 1+__kernel(y,n,n);

    /* assemble symmetric matrix */
    // A->assemble_Matrix(diagonal, tolerance, 's');

    /* compute factorization */
    // A->compute_Factor();

    /* compute determinant */
    // A->compute_Determinant(logdet);

    // Computing log-determinant using Cholesky:
    Eigen::LLT<MatrixXd> P;
    P.compute(B);
    for(int i=0; i<P.matrixL().rows(); ++i)
    {
        logdet += log(P.matrixL()(i,i));
    }
    logdet = 2 * logdet;

    printf("exact: %.14g\n", exact);
    printf("HODLR: %.14g\n", logdet);
    printf("relative error: %g\n", 1-fabs(logdet/exact));

    return 0;
}
