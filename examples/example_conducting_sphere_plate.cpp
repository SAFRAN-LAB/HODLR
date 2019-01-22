#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

using std::setprecision;

// This kernel arises in the Casimir effect between a sphere of radius R and a
// plate separated by a distance L. The kernel corresponds to the
// high-temperature limit for Drude metals and m=0, see also Ref 1.
//
// References:
// [1] Giuseppe Bimonte, Classical Casimir interaction of a perfectly
//     conducting sphere and plate, Phys. Rev. D 95, 065004 (2017)
//     https://doi.org/10.1103/PhysRevD.95.065004
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

    // The matrix is given by:
    //       D = Id - M_ij
    //  where
    //       M_ij = y^(i+j+1) * (i+j)!/(i! j!)
    //  with the indices 1 <= i,j <= ldim, and
    //       y = 0.5*R/(R+L).
    double getMatrixEntry(int i, int j) 
    {
        const int l1 = i+1, l2 = j+1;
        const double delta = (i == j) ? 1 : 0; // Kronecker delta

        return delta-exp((l1 + l2 + 1) * y + lgamma(1 + l1 + l2) - lgamma(1 + l1) - lgamma(1 + l2));
    };
};

int main(int argc, char *argv[]) 
{
    // exact (analytical) result from Ref. [1]
    const double exact = -117.00871690447;

    // RbyL: aspect ratio; R corresponds to the radius of the sphere and L to
    // the separation between sphere and plane
    const double RbyL = 42000;

    // The matrices have infinite dimension; truncate matrix according to
    // dim = 10*R/L
    const unsigned int dim = 10 * RbyL;

    // tolerance
    const double tolerance = 1e-11;

    // nLeaf is the size (number of rows of the matrix) of the smallest block
    // at the leaf level. The number of levels in the tree is given by
    // n_levels=log_2(N/nLeaf) where N denotes the dimension of the matrix.
    const int nLeaf = 100;
    int n_levels  = log(dim / nLeaf) / log(2);

    kernel K(dim, RbyL);

    HODLR_Tree* T = new HODLR_Tree(n_levels, tolerance, &K);

    std::cout << "========================= Problem Parameters =========================" << std::endl;
    std::cout << "Matrix Size                        :" << dim << std::endl;
    std::cout << "Leaf Size                          :" << nLeaf << std::endl;
    std::cout << "Number of Levels in Tree           :" << n_levels << std::endl;
    std::cout << "Tolerance                          :" << tolerance << std::endl << std::endl;

    // Assemble symmetric matrix 
    // If we are assembling a symmetric matrix:
    bool is_sym = true;
    // If we know that the matrix is also PD:
    // By toggling this flag to true, the factorizations are performed using Cholesky
    // Useful when you want the factorization as WW^T 
    bool is_pd = true;

    T->assembleTree(is_sym, is_pd);
    // T->plotTree("plot.svg");

    // Compute factorization 
    T->factorize();

    // Compute determinant 
    double logdet = T->logDeterminant();

    std::cout << "Log determinant(Exact):" << setprecision(16) << exact << std::endl;
    std::cout << "Log determinant(HODLR):" << setprecision(16) << logdet << std::endl;
    std::cout << "Relative error is: "  << setprecision(16) << fabs(1 - logdet / exact) << std::endl;

    return 0;
}
