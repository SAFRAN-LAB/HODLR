#include "HODLR_Matrix.hpp"

class Kernel : public HODLR_Matrix
{
private:
    Vec x, y;

public:
    // Constructor:
    Kernel(int N);
    // Destructor:
    ~Kernel();

    // Returns a particular entry of the matrix:
    dtype getMatrixEntry(int i, int j) const;

    // Returns Matrix x:
    // dtype
};
