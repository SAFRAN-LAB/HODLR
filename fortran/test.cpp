#include "Kernel.hpp"

int main(int argc, char* argv[])
{
    Kernel* K = new Kernel(10);
    std::cout << K->getMatrixEntry(0, 0) << std::endl;
}
