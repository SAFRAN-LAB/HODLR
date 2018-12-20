#include "HODLR_Tree.hpp"

void HODLR_Tree::factorizeLeafSPD(int k) 
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
        r       = tree[l][parent]->rank[child]; // NOTE: Here the rank for both childs is the same

        tree[l][parent]->Q_factor[child].block(t_start, 0, size, r) =   
        this->solveLeafSymmetricFactor(k, tree[l][parent]->Q_factor[child].block(t_start, 0, size, r));
    }
}

void HODLR_Tree::qr(int j, int k)
{
    int min0 = std::min(tree[j][k]->Q_factor[0].rows(), 
                        tree[j][k]->Q_factor[0].cols()
                       );
    
    int min1 = std::min(tree[j][k]->Q_factor[1].rows(),
                        tree[j][k]->Q_factor[1].cols()
                       );

    Eigen::HouseholderQR<Eigen::MatrixXd> qr(tree[j][k]->Q_factor[0]);

    // Getting thin Q:
    tree[j][k]->Q_factor[0] = qr.householderQ() * 
                              Eigen::MatrixXd::Identity(tree[j][k]->Q_factor[0].rows(), min0);
    // K = I * R
    tree[j][k]->K          *= qr.matrixQR().block(0, 0, min0, 
                                                  tree[j][k]->Q_factor[0].cols()
                                                 ).triangularView<Eigen::Upper>();

    qr = Eigen::HouseholderQR<Eigen::MatrixXd>(tree[j][k]->Q_factor[1]);
    // Getting thin Q:
    tree[j][k]->Q_factor[1] = qr.householderQ() * 
                              Eigen::MatrixXd::Identity(tree[j][k]->Q_factor[1].rows(), min1);
    // K = I * R * RT
    tree[j][k]->K          *= qr.matrixQR().block(0, 0, min1, 
                                                  tree[j][k]->Q_factor[1].cols()
                                                 ).triangularView<Eigen::Upper>().transpose();
}

// Obtains the QR decomposition of all nodes at this level:
void HODLR_Tree::qrForLevel(int level)
{
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[level]; k++)
    {
        qr(level, k);
    }
}

void HODLR_Tree::factorizeNonLeafSPD(int j, int k) 
{
    int r0 = tree[j][k]->rank[0];
    int r1 = tree[j][k]->rank[1];

    tree[j][k]->K_factor_LLT.compute(  MatrixXd::Identity(r0, r1) 
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

        if(tree[l][parent]->Q_factor[child].cols() > 0)
        {
            tree[l][parent]->Q_factor[child].block(t_start, 0, size, r) =   
            this->solveNonLeafSymmetricFactor(j, k, tree[l][parent]->Q_factor[child].block(t_start, 0, size, r));
        }
    }
}

void HODLR_Tree::factorizeSPD() 
{
    // Initializing for the non-leaf levels:
    for(int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            int &r0 = tree[j][k]->rank[0];
            int &r1 = tree[j][k]->rank[1];

            // Initializing the factorized matrices for the left and right child:
            tree[j][k]->Q_factor[0] = tree[j][k]->Q[0];
            tree[j][k]->Q_factor[1] = tree[j][k]->Q[1];
            tree[j][k]->K           = MatrixXd::Identity(r0, r1); // NOTE: r0 == r1
        }
    }

    // Factorizing the leaf levels:
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        this->factorizeLeafSPD(k);
    }

    qrForLevel(n_levels - 1);
    // Factorizing the nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            this->factorizeNonLeafSPD(j, k);
        }

        qrForLevel(j-1);
    }

    if(n_levels > 0)
    {
        tree[0][0]->K_factor_LLT.compute(  MatrixXd::Identity(tree[0][0]->rank[0], tree[0][0]->rank[1]) 
                                         - tree[0][0]->K.transpose() * tree[0][0]->K
                                        );
    }
}

MatrixXd HODLR_Tree::solveLeafSymmetricFactor(int k, MatrixXd b) 
{
    MatrixXd x = tree[n_levels][k]->K_factor_LLT.matrixL().solve(b);
    return x;
}

MatrixXd HODLR_Tree::solveLeafSymmetricFactorTranspose(int k, MatrixXd b) 
{
    MatrixXd x = tree[n_levels][k]->K_factor_LLT.matrixL().transpose().solve(b);
    return x;
}

MatrixXd HODLR_Tree::solveNonLeafSymmetricFactor(int j, int k, MatrixXd b) 
{
    int n0 = tree[j][k]->Q[0].rows();
    int n1 = tree[j][k]->Q[1].rows();
    int r  = b.cols();

    MatrixXd tmp = tree[j][k]->Q_factor[1].transpose() * b.block(n0, 0, n1, r);
    // TODO:Not clear about this yet
    b.block(n0, 0, n1, r) -=   tree[j][k]->Q_factor[1]
                                 * (  tree[j][k]->K_factor_LLT.matrixL().solve((tree[j][k]->K.transpose()
                                    *(tree[j][k]->Q_factor[0].transpose() * b.block(0, 0, n0, r))) - tmp) + tmp
                                   );
    return(b);
}

MatrixXd HODLR_Tree::solveNonLeafSymmetricFactorTranspose(int j, int k, Eigen::MatrixXd b)
{
    int n0 = tree[j][k]->Q[0].rows();
    int n1 = tree[j][k]->Q[1].rows();
    int r  = b.cols();

    // TODO:Not clear about this yet
    Eigen::MatrixXd xtmp = tree[j][k]->Q_factor[1].transpose() * b.block(n0,0,n1,b.cols());
    Eigen::MatrixXd ytmp = tree[j][k]->K_factor_LLT.matrixL().transpose().solve(xtmp);
    
    b.block(0 , 0, n0, r) -= tree[j][k]->Q_factor[0] * (tree[j][k]->K * ytmp);
    b.block(n0, 0, n1, r) -= tree[j][k]->Q_factor[1] * (xtmp - ytmp);
    
    return(b);
}

MatrixXd HODLR_Tree::solveSymmetricFactor(MatrixXd b)
{
    int start, size;
    MatrixXd x = MatrixXd::Zero(b.rows(),b.cols());
    
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
    }

    return x;
} 

MatrixXd HODLR_Tree::solveSymmetricFactorTranspose(MatrixXd b) 
{
    int start, size;
    MatrixXd x = MatrixXd::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    // Factoring out the leaf nodes:
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = 
        this->solveLeafSymmetricFactorTranspose(k, b.block(start, 0, size, r));
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
            this->solveNonLeafSymmetricFactorTranspose(j, k, b.block(start, 0, size, r));
        }
    }

    return x;
}

MatrixXd HODLR_Tree::solveSPD(MatrixXd b) 
{   
    return(solveSymmetricFactorTranspose(solveSymmetricFactor(b)));
}

double HODLR_Tree::logDeterminantSPD()
{
    double log_det = 0.0;

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
