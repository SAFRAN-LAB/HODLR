import numpy as np
import pyhodlrlib

# Size of the matrix:
N = 500
x = np.sort(np.random.rand(N))

# Returning the Gaussian Kernel:
class Kernel(pyhodlrlib.HODLR_Matrix):
    def getMatrixEntry(self, i, j):
        if(i == j):
            return 10
        else:
            return np.exp(-(x[i] - x[j])**2)

def test_functions(factorization_type, is_sym, is_pd):
    K = Kernel(N)
    # What we are doing here is explicitly generating 
    # the matrix from its entries
    A = K.getMatrix(0, 0, N, N)
    # Tolerance for all factorizations:
    eps = 1e-12
    # Declaring the Factorizer Object:
    F = pyhodlrlib.Matrix_Factorizer(K, factorization_type)

    # Size of leaf level:
    M = 50
    # Number of levels:
    n_levels = int(np.log(N / M) / np.log(2))

    # Creating the HODLR Tree object:
    T = pyhodlrlib.HODLR_Tree(n_levels, eps, F)
    T.assembleTree(is_sym, is_pd)

    # Random vector to take product with:
    x = np.random.rand(N)
    # Finding b using HODLR:
    b_hodlr = T.matmatProduct(x)
    # Finding b using direct MatVec:
    b = A @ x
    # Verifying the accuracy of the MatMat:
    assert(np.linalg.norm(b_hodlr.ravel() - b.ravel()) / np.linalg.norm(b) < N * eps)

    # Factorize elements of the tree:
    T.factorize()
    # Solving for x in A x = b:
    x_hodlr = T.solve(b)
    # Computing the relative error:
    assert(np.linalg.norm(x_hodlr.ravel() - x) /  np.linalg.norm(x) < N * eps)

    # Finding log determinant:
    logdet_hodlr = T.logDeterminant()
    # Finding logdet:
    logdet = 2 * np.sum(np.log(abs(np.diag(np.linalg.cholesky(A)))))
    assert(abs(logdet_hodlr - logdet) / logdet < N * eps)

    # When system described is SPD:
    if(is_sym and is_pd):
        # Setting y = W^T x
        y = T.symmetricFactorTransposeProduct(x)
        # Getting b = W (W^T x)
        b_hodlr = T.symmetricFactorProduct(y)
        assert(np.linalg.norm(b_hodlr.ravel() - b.ravel()) / np.linalg.norm(b) <  N * eps)

        # Directly obtaining the symmetric factor matrix:
        W = T.getSymmetricFactor()
        assert(np.mean(abs(W @ W.T - A)) < N * eps)

    return

test_functions('rookPivoting', False, False)
test_functions('rookPivoting', True, False)
test_functions('rookPivoting', True, True)

test_functions('queenPivoting', False, False)
test_functions('queenPivoting', True, False)
test_functions('queenPivoting', True, True)

test_functions('SVD', False, False)
test_functions('SVD', True, False)
test_functions('SVD', True, True)

test_functions('RRQR', False, False)
test_functions('RRQR', True, False)
test_functions('RRQR', True, True)

print('Reached End of Test File Successfully! All functions work as intended!')
