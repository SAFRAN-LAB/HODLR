import numpy as np
import pyhodlrlib

# Size of the matrix:
N = 100

x = 3 + np.random.rand(N)
y = 7 + np.random.rand(N)

# Returning the Gaussian Kernel:
class Kernel(pyhodlrlib.HODLR_Matrix):
    def getMatrixEntry(self, i, j):
        # Gaussian Kernel:
        return np.exp(-(x[i] - y[j])**2)

K = Kernel(N)
A = K.getMatrix(0, 0, N, N)
# print('Matrix A is:\n', A)

# Tolerance for all factorizations:
eps = 1e-16

# Declaring the Factorizer Object:
F = pyhodlrlib.Matrix_Factorizer(K, 'rookPivoting')
# Getting the factorization of Matrix A directly:
F.factorize(1e-12, 0, 0, N, N)
L = F.getL()
R = F.getR()
# Verifying the accuracy of the factorization:
print('Error in Factorization:', np.sum(abs(A - L @ R.T)))

# Size of leaf level:
M = 100

# Number of levels:
n_levels = int(np.log(N / M) / np.log(2))

# Creating the HODLR Tree object:
T = pyhodlrlib.HODLR_Tree(n_levels, eps, F)
T.assembleTree(False, False)

# Random vector to take product with:
x = np.random.rand(N, N)

# Finding b using HODLR:
b_hodlr = T.matmatProduct(x)
# Finding b using direct MatMat:
b = A @ x
# Verifying the accuracy of the MatMat:
print('Error in Matrix-Matrix Multiplication:', np.sum(abs(b_hodlr - b)))

# Factorize elements of the tree:
T.factorize()

# Solving for x in A x = b:
x_hodlr = T.solve(b)
# Verifying the accuracy of solve:
print('Error in Solve:', np.sum(abs(x_hodlr - x)))

# Setting y = W^T x
y = T.symmetricFactorTransposeProduct(x)
# Getting b = W (W^T x)
b_hodlr = T.symmetricFactorProduct(y)
print('Error in Symmetric Factor Multiplications:', np.sum(abs(b_hodlr - b)))

# Finding log determinant:
logdet_hodlr = T.logDeterminant()
print('Error in Determinant Computation:', np.sum(abs(logdet_hodlr - np.log(np.linalg.det(A)))))

# Directly obtaining the symmetric factor matrix:
W = T.getSymmetricFactor()
print('Error in Getting Symmetric Factor:', np.sum(abs(W @ W.T - A)))
