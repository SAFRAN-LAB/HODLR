#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

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
            return 10;
        }
        
        // Otherwise:
        else
        {   
            dtype R2;
            // Initializing:
            R2 = 0;

            for(int k = 0; k < dim; k++)
            {
                R2 += (x(i,k) - x(j,k)) * (x(i,k) - x(j,k));
            }

            // Gaussian Kernel: e(-R^2)
            return exp(-R2);
        }
    }

    // Destructor:
    ~Kernel(){}
};
