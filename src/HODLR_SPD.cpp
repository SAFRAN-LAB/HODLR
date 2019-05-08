#include "HODLR.hpp"

// Factorizing out the leaf level matrices:
// This is the first step of the process of factoring out:
// That is we are making K = K_κ * K^(κ-1) * K_κ^T
void HODLR::factorizeLeafSPD(int k) 
{
    int child;
    int parent = k;
    int size   = tree[n_levels][k]->n_size;
    
    int t_start, r;
    
    tree[n_levels][k]->K_factor_LLT.compute(tree[n_levels][k]->K);
    for(int l = n_levels - 1; l >= 0; l--) 
    {
        child   = parent % 2;
        parent  = parent / 2;
        t_start = tree[n_levels][k]->n_start - tree[l][parent]->c_start[child];
        r       = tree[l][parent]->rank[child]; // NOTE: Here the rank for both children is the same

        // We factor out the leaf level by applying inv(L) to the appropriate subblock:
        tree[l][parent]->Q[child].block(t_start, 0, size, r) =   
        this->solveLeafSymmetricFactor(k, tree[l][parent]->Q[child].block(t_start, 0, size, r));
    }
}

// Applies QR factorization of the node at the given level and node number:
// Logic here is that, we have express U0 * V1T instead as:
// U0 * (V1)^T = (Q0 * R0) * (Q1 * R1)^T = Q0 * R0 * R1^T * Q1^T 
// = Q0 * K * Q1^T; where K = R0 * R1^T 
void HODLR::qr(int j, int k)
{   
    // We are finding the minimum size along the dimension since we want to
    // get the thinQ of the QR factorization for U1 and V2:
    
    // min0 = min(N, r)
    int min0 = std::min(tree[j][k]->Q[0].rows(), 
                        tree[j][k]->Q[0].cols()
                       );
    
    // min0 = min(N, r)
    int min1 = std::min(tree[j][k]->Q[1].rows(),
                        tree[j][k]->Q[1].cols()
                       );

    // Performing QR factorization of U0, U0 = Q0 * R0:
    Eigen::HouseholderQR<Mat> qr(tree[j][k]->Q[0]);

    // Getting thin Q:
    // householderQ has shape (N, N)
    // multiplying we can make that (N, r)
    tree[j][k]->Q[0] = qr.householderQ() * 
                       Mat::Identity(tree[j][k]->Q[0].rows(), min0);

    // K = R0
    tree[j][k]->K = qr.matrixQR().block(0, 0, min0, min0).triangularView<Eigen::Upper>();

    // Performing QR factorization: V1 = Q1 * R1
    qr = Eigen::HouseholderQR<Mat>(tree[j][k]->Q[1]);
    // Getting thin Q:
    tree[j][k]->Q[1] = qr.householderQ() * 
                       Mat::Identity(tree[j][k]->Q[1].rows(), min1);
    // K = R0 * R1T
    tree[j][k]->K *= qr.matrixQR().block(0, 0, min1, min1).triangularView<Eigen::Upper>().transpose();
}

// Obtains the QR decomposition of all nodes at this level:
void HODLR::qrForLevel(int level)
{
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[level]; k++)
    {
        qr(level, k);
    }
}

void HODLR::factorizeNonLeafSPD(int j, int k) 
{
    int r0 = tree[j][k]->rank[0];
    int r1 = tree[j][k]->rank[1];

    // Computes L * L^T = (I - K^T * K)
    //                  = (I - (R0 * R1^T)^T * (R0 * R1^T))
    //                  = (I - (R1 * R0^T)   * (R0 * R1^T))
    tree[j][k]->K_factor_LLT.compute(  Mat::Identity(r0, r1) 
                                     - tree[j][k]->K.transpose() * tree[j][k]->K
                                    );

    int parent = k;
    int child  = k;
    int size   = tree[j][k]->n_size;
    int t_start, r;

    for(int l = j - 1; l >= 0; l--) 
    {
        child   = parent % 2;
        parent  = parent / 2;
        t_start = tree[j][k]->n_start - tree[l][parent]->c_start[child];
        r       = tree[l][parent]->rank[child];

        if(tree[l][parent]->Q[child].cols() > 0)
        {
            tree[l][parent]->Q[child].block(t_start, 0, size, r) =   
            this->solveNonLeafSymmetricFactor(j, k, tree[l][parent]->Q[child].block(t_start, 0, size, r));
        }
    }
}

void HODLR::factorizeSPD()
{
    // Initializing for the non-leaf levels:
    for(int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            // Initializing the factorized matrices for the left and right child:
            // In the symmetric case, the HODLR matrix is given by:
            // |           |           |
            // |    K_0    |  U0 * V1T |
            // |           |           |
            // -------------------------
            // |           |           |
            // |  V1 * U0T |    K_1    |
            // |           |           |
            // -------------------------
            // We initialize Qs to be U0 and V1. In the following steps, 
            // we will get the QR factorization of these matrices:
            tree[j][k]->Q[0] = tree[j][k]->U[0];
            tree[j][k]->Q[1] = tree[j][k]->V[1];
        }
    }

    // Factorizing the leaf levels:
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        this->factorizeLeafSPD(k);
    }

    // Factorizing the nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        qrForLevel(j);
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            this->factorizeNonLeafSPD(j, k);
        }
    }
}

// Returns inv(L) * b 
Mat HODLR::solveLeafSymmetricFactor(int k, Mat b) 
{
    Mat x = tree[n_levels][k]->K_factor_LLT.matrixL().solve(b);
    return x;
}

// Returns inv(LT) * b 
Mat HODLR::solveLeafSymmetricFactorTranspose(int k, Mat b) 
{
    Mat x = tree[n_levels][k]->K_factor_LLT.matrixL().transpose().solve(b);
    return x;
}

Mat HODLR::solveNonLeafSymmetricFactor(int j, int k, Mat b) 
{
    int n0 = tree[j][k]->U[0].rows();
    int n1 = tree[j][k]->V[1].rows();
    int r  = b.cols();

    // tmp = Q1^T * b
    Mat tmp = tree[j][k]->Q[1].transpose() * b.block(n0, 0, n1, r);
    // What we are trying to solve is:
    // (I + Q1 * K^T * Q0^T) * x = b
    // by Sherman Morrisson Woodbury Formula:
    // x = b - Q1 * inv(inv(K^T) + Q0^T * Q1) * Q1^T
    b.block(n0, 0, n1, r) -=   tree[j][k]->Q[1]
                             * (  tree[j][k]->K_factor_LLT.matrixL().solve((tree[j][k]->K.transpose()
                                *(tree[j][k]->Q[0].transpose() * b.block(0, 0, n0, r))) - tmp) + tmp
                               );
    return(b);
}

Mat HODLR::solveNonLeafSymmetricFactorTranspose(int j, int k, Mat b)
{
    int n0 = tree[j][k]->U[0].rows();
    int n1 = tree[j][k]->V[1].rows();
    int r  = b.cols();

    // xtmp = Q1T * b
    // ytmp = inv(LT) * Q1T * b
    Mat xtmp = tree[j][k]->Q[1].transpose() * b.block(n0, 0, n1, r);
    Mat ytmp = tree[j][k]->K_factor_LLT.matrixL().transpose().solve(xtmp);
    
    // b = b - Q0 * R * RT * inv(L)
    // b = b - Q1 * (ytmp - ytmp)
    b.block(0 , 0, n0, r) -= tree[j][k]->Q[0] * (tree[j][k]->K * ytmp);
    b.block(n0, 0, n1, r) -= tree[j][k]->Q[1] * (xtmp - ytmp);
    
    return(b);
}

Mat HODLR::solveSymmetricFactor(Mat b)
{
    int start, size;
    Mat x = Mat::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    // Factoring out the leaf nodes:
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = 
        this->solveLeafSymmetricFactor(k, b.block(start, 0, size, r));
    }

    b = x;
    
    // Factoring out over nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            start = tree[j][k]->n_start;
            size  = tree[j][k]->n_size;
            
            x.block(start, 0, size, r) = 
            this->solveNonLeafSymmetricFactor(j, k, b.block(start, 0, size, r));
        }

        b = x;
    }

    return x;
} 

Mat HODLR::solveSymmetricFactorTranspose(Mat b) 
{
    int start, size;
    Mat x = Mat::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    // Factoring out over nonleaf levels:
    for(int j = 0; j < n_levels; j++) 
    {
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            start = tree[j][k]->n_start;
            size  = tree[j][k]->n_size;
            
            x.block(start, 0, size, r) = 
            this->solveNonLeafSymmetricFactorTranspose(j, k, b.block(start, 0, size, r));
        }

        b = x;
    }

    // Factoring out the leaf nodes:
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = 
        this->solveLeafSymmetricFactorTranspose(k, b.block(start, 0, size, r));
    }

    return x;
}

Mat HODLR::solveSPD(Mat b) 
{   
    return(solveSymmetricFactorTranspose(solveSymmetricFactor(b)));
}

Mat HODLR::SymmetricFactorNonLeafProduct(int j, int k, Mat b) 
{
    int n0                        = tree[j][k]->U[0].rows();
    int n1                        = tree[j][k]->V[1].rows();
    Mat tmp                       = tree[j][k]->Q[1].transpose() * b.block(n0, 0, n1, b.cols());
    b.block(n0, 0, n1, b.cols()) += tree[j][k]->Q[1]*(  (  tree[j][k]->K.transpose() 
                                                         * tree[j][k]->Q[0].transpose() 
                                                         * b.block(0, 0, n0, b.cols())
                                                        ) 
                                                      + (  (Mat)tree[j][k]->K_factor_LLT.matrixL() 
                                                         - Mat::Identity(tree[j][k]->rank[0], tree[j][k]->rank[1])
                                                        ) * tmp
                                                     );

    return(b);
}

Mat HODLR::SymmetricFactorTransposeNonLeafProduct(int j, int k, Mat b)
{
    int n0                        = tree[j][k]->U[0].rows();
    int n1                        = tree[j][k]->V[1].rows();
    Mat tmp                       = tree[j][k]->Q[1].transpose() * b.block(n0, 0, n1, b.cols());
    b.block(0,  0, n0, b.cols()) += tree[j][k]->Q[0] * tree[j][k]->K * tmp;
    b.block(n0, 0, n1, b.cols()) += tree[j][k]->Q[1] * ((  (Mat)tree[j][k]->K_factor_LLT.matrixL().transpose() 
                                                          - Mat::Identity(tree[j][k]->rank[0], tree[j][k]->rank[1])
                                                        ) * tmp
                                                       );

    return(b);
}

Mat HODLR::symmetricFactorProduct(Mat b)
{
    int start, size;
    Mat x = Mat::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    for(int j = 0; j < n_levels; j++)
    {
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            start = tree[j][k]->n_start;
            size  = tree[j][k]->n_size;
            
            x.block(start, 0, size, r) = 
            this->SymmetricFactorNonLeafProduct(j, k, b.block(start, 0, size, r));
        }

        b = x;
    } 
    
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = 
        tree[n_levels][k]->K_factor_LLT.matrixL() * b.block(start, 0, size, r);
    }

    return x;
}

Mat HODLR::symmetricFactorTransposeProduct(Mat b)
{
    int start, size;
    Mat x = Mat::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = 
        tree[n_levels][k]->K_factor_LLT.matrixL().transpose() * b.block(start, 0, size, r);
    }

    b = x;
    
    // Factoring out over nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            start = tree[j][k]->n_start;
            size  = tree[j][k]->n_size;
            
            x.block(start, 0, size, r) = 
            this->SymmetricFactorTransposeNonLeafProduct(j, k, b.block(start, 0, size, r));
        }

        b = x;
    }

    return x;
}

Mat HODLR::getSymmetricFactor()
{
    if(n_levels == 0)
        return tree[0][0]->K_factor_LLT.matrixL();

    Mat Wc;
    Mat Rc = Mat::Identity(N,N);

    for(int j = 0; j < n_levels; j++)
    {
        Wc = Mat::Zero(N, N);
        
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++)
        {
            int r      = tree[j][k]->rank[0];
            int n0     = tree[j][k]->Q[0].rows();
            int n1     = tree[j][k]->Q[1].rows();
            Mat T = Mat::Identity(n0 + n1, n0 + n1);
            
            T.block(n0, 0, n1, n0)   = tree[j][k]->Q[1] * (  tree[j][k]->K.transpose()
                                                           * tree[j][k]->Q[0].transpose()
                                                          );

            T.block(n0, n0, n1, n1) += tree[j][k]->Q[1]*(  (  (Mat)tree[j][k]->K_factor_LLT.matrixL() 
                                                            -  Mat::Identity(r, r)
                                                           ) 
                                                         * tree[j][k]->Q[1].transpose()
                                                        );

            Wc.block(tree[j][k]->n_start, tree[j][k]->n_start, tree[j][k]->n_size, tree[j][k]->n_size) = T;
        }

        Rc = Wc * Rc;
    }

    Wc = Mat::Zero(N, N);

    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++)
    {
        Wc.block(tree[n_levels - 1][k / 2]->c_start[k % 2], 
                 tree[n_levels - 1][k / 2]->c_start[k % 2], 
                 tree[n_levels - 1][k / 2]->c_size[k % 2], 
                 tree[n_levels - 1][k / 2]->c_size[k % 2]) = tree[n_levels][k]->K_factor_LLT.matrixL();
    }

    Rc = Wc * Rc;
    return Rc;
}

dtype HODLR::logDeterminantSPD()
{
    dtype log_det = 0.0;

    for(int j = n_levels; j >= 0; j--) 
    {
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            if(tree[j][k]->K.size() > 0)
            {
                for(int l = 0; l < tree[j][k]->K_factor_LLT.matrixL().rows(); l++) 
                {   
                    log_det += log(fabs(tree[j][k]->K_factor_LLT.matrixL()(l,l)));
                }
            }
        }
    }

    log_det *= 2;
    return(log_det);
}
