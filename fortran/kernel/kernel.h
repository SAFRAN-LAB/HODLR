#ifdef __cplusplus // Are we compiling this with a C++ compiler ?
extern "C" 
{
    class Kernel;
    typedef Kernel KERNEL;
#else
    // From the C side, we use an opaque pointer.
    typedef struct KERNEL KERNEL;
#endif

// Constructor
KERNEL* create_kernel(int N);

// Destructor
void delete_kernel(KERNEL* K);

// The const qualificators maps from the member function 
// to pointers to the class instances.
double get_matrix_entry(const KERNEL* K, int i, int j);

double get_vector_x(const KERNEL* K);

#ifdef __cplusplus
}
#endif
