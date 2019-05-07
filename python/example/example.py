import numpy as np
import pyhodlrlib

# Size of the matrix:
N = 1000
x = np.sort(np.random.rand(N))
# Size of leaf level:
M = 200

# Returning the Gaussian Kernel:
class Kernel(pyhodlrlib.HODLR_Matrix):
    def getMatrixEntry(self, i, j):
        if(i == j):
            return 10
        else:
            return np.exp(-(x[i] - x[j])**2)

K = Kernel(N)
# What we are doing here is explicitly generating 
# the matrix from its entries
A = K.getMatrix(0, 0, N, N)
# Tolerance for all factorizations:
eps = 1e-12

# If we are assembling a symmetric matrix:
is_sym = True
# If we know that the matrix is also PD:
# By setting the matrix to be symmetric-positive definite, 
# we trigger the fast symmetric factorization method to be used
# In all other cases the fast factorization method is used
is_pd = True
# Creating the HODLR object:
T = pyhodlrlib.HODLR(N, M, eps, K, 'rookPivoting', is_sym, is_pd)

# Random vector to take product with:
x = np.random.rand(N)
# Finding b using HODLR:
b_hodlr = T.matmatProduct(x)
# Finding b using direct MatVec:
b = A @ x
# Verifying the accuracy of the MatMat:
print('Error in Matrix-Matrix Multiplication:', np.linalg.norm(b_hodlr.ravel() - b.ravel()) / np.linalg.norm(b))

# Factorize elements of the tree:
T.factorize()
# Solving for x in A x = b:
x_hodlr = T.solve(b)
# Computing the relative error:
print('Error in Solve:', np.linalg.norm(x_hodlr.ravel() - x) /  np.linalg.norm(x))

# Finding log determinant:
logdet_hodlr = T.logDeterminant()
# Finding logdet:
logdet = 2 * np.sum(np.log(abs(np.diag(np.linalg.cholesky(A)))))
print('Error in Log-Determinant Computation:', abs(logdet_hodlr - logdet))

# When system described is SPD:
if(is_sym and is_pd):
    # Setting y = W^T x
    y = T.symmetricFactorTransposeProduct(x)
    # Getting b = W (W^T x)
    b_hodlr = T.symmetricFactorProduct(y)
    print('Error in Symmetric Factor Multiplications:', np.linalg.norm(b_hodlr.ravel() - b.ravel()) / np.linalg.norm(b))

    # Directly obtaining the symmetric factor matrix:
    W = T.getSymmetricFactor()
    print('Error in Getting Symmetric Factor:', np.mean(abs(W @ W.T - A)))
