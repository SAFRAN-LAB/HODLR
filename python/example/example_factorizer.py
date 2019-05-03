import numpy as np
import pyhodlrlib

# Size of the matrix:
N = 1000
x = 3 + np.sort(np.random.rand(N))
y = 7 + np.sort(np.random.rand(N))

# Returning the Gaussian Kernel:
class Kernel(pyhodlrlib.HODLR_Matrix):
    def getMatrixEntry(self, i, j):
        return np.exp(-(x[i] - y[j])**2)

K = Kernel(N)
# What we are doing here is explicitly generating 
# the matrix from its entries
A = K.getMatrix(0, 0, N, N)
# Tolerance for all factorizations:
eps = 1e-12

# Declaring the Factorizer Object:
F = pyhodlrlib.Matrix_Factorizer(K, 'rookPivoting')
# Getting the factorization of Matrix A directly:
F.factorize(eps, 0, 0, N, N)
L = F.getL()
R = F.getR()
# Verifying the accuracy of the factorization:
print('Error in Factorization:', np.mean(abs(A - L @ R.T)))
