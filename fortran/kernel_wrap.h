// Are we compiling this with a C++ compiler?
#ifdef __cplusplus 
extern "C" 
{
    class Kernel;
    typedef Kernel KERNEL;
#else
    // From the C side, we use an opaque pointer.
    typedef struct KERNEL KERNEL;
#endif

// Constructor
KERNEL* createKernel(int N);

// Destructor
void deleteKernel(KERNEL* K);

// The const qualificators maps from the member function to pointers to the instance:
double* getMatrix(const KERNEL* K, int j, int k, int n_rows, int n_cols);

#ifdef __cplusplus
}
#endif
