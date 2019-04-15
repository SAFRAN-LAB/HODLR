#include "HODLR_Matrix.hpp"

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
    };
    
    dtype getMatrixEntry(int i, int j) 
    {
        return 1 / abs((x(i)-y(j)));
    }

    // Destructor:
    ~Kernel() 
    {
        std::cout << "Kernel Object Destroyed!!" << std::endl;
    };
};
