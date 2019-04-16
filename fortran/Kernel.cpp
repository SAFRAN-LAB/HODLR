#include "Kernel.hpp"

Kernel::Kernel(int N) : HODLR_Matrix(N)
{
    srand(time(NULL));
    x = 2 * Vec::Ones(N) + Vec::Random(N);
    y = 7 * Vec::Ones(N) + Vec::Random(N);

    std::sort(x.data(), x.data() + x.size(), std::greater<double>());
    std::sort(y.data(), y.data() + y.size());
}

Kernel::~Kernel()
{
    std::cout << "Kernel Object Destroyed!!" << std::endl;
}

dtype Kernel::getMatrixEntry(int i, int j) const
{
    return 1 / abs((x(i)-y(j)));
}
