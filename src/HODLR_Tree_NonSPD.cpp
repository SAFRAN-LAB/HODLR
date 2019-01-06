#include "HODLR_Tree.hpp"

// Factorizing out the leaf nodes:
// This is the first step of the process of factoring out:
// That is we are making K^(κ) = K_κ * K^(κ-1)
void HODLR_Tree::factorizeLeafNonSPD(int k) 
{
    int child;
    int parent = k;
    int size   = tree[n_levels][k]->n_size;
    
    int t_start, r;
    
    tree[n_levels][k]->K_factor_LU.compute(tree[n_levels][k]->K);
    for(int l = n_levels - 1; l >= 0; l--) 
    {
        child   = parent % 2;
        parent  = parent / 2;
        t_start = tree[n_levels][k]->n_start - tree[l][parent]->c_start[child];
        r       = tree[l][parent]->rank[child];

        // We factor out the leaf level by applying inv(K) to the appropriate subblock:
        tree[l][parent]->U_factor[child].block(t_start, 0, size, r) =   
        this->solveLeafNonSPD(k, tree[l][parent]->U_factor[child].block(t_start, 0, size, r));
    }
}

// Factorizing out the nonleaf nodes:
// That is we are making K^(k) = K_k * K^(k-1)
void HODLR_Tree::factorizeNonLeafNonSPD(int j, int k) 
{
    int r0 = tree[j][k]->rank[0];
    int r1 = tree[j][k]->rank[1];

    if(r0 > 0 || r1 > 0)
    {
        // Through the steps below, we are making:
        //     |           |           |
        //     |     I     |  V0T * U0 |
        //     |           |           |
        // K = ------------------------- = I + V^T U
        //     |           |           |
        //     |  V1T * U1 |     I     |
        //     |           |           |
        //     -------------------------

        tree[j][k]->K.block(0, r0, r0, r1)  =   
        tree[j][k]->V_factor[1].transpose() * tree[j][k]->U_factor[1];

        tree[j][k]->K.block(r0, 0, r1, r0)  =   
        tree[j][k]->V_factor[0].transpose() * tree[j][k]->U_factor[0];

        // Computing LU factorization of K:
        tree[j][k]->K_factor_LU.compute(tree[j][k]->K);

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

            // This updates U's by using the Sherman Morrisson Woodbury formula:
            if(tree[l][parent]->U_factor[child].cols() > 0)
            {
                tree[l][parent]->U_factor[child].block(t_start, 0, size, r) =   
                this->solveNonLeafNonSPD(j, k, tree[l][parent]->U_factor[child].block(t_start, 0, size, r));
            }
        }
    }
}

void HODLR_Tree::factorizeNonSPD() 
{
    // Initializing for the non-leaf levels:
    for(int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            int &r0 = tree[j][k]->rank[0];
            int &r1 = tree[j][k]->rank[1];

            tree[j][k]->U_factor[0] = tree[j][k]->U[0];
            tree[j][k]->U_factor[1] = tree[j][k]->U[1];
            tree[j][k]->V_factor[0] = tree[j][k]->V[0];
            tree[j][k]->V_factor[1] = tree[j][k]->V[1];
            tree[j][k]->K           = Mat::Identity(r0 + r1, r0 + r1);
        }
    }

    // Factorizing the leaf levels:
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        this->factorizeLeafNonSPD(k);
    }

    // Factorizing the nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            this->factorizeNonLeafNonSPD(j, k);
        }
    }
}

// Solve at the leaf is just directly performed by solving Kx = b:
Mat HODLR_Tree::solveLeafNonSPD(int k, Mat b) 
{
    Mat x = tree[n_levels][k]->K_factor_LU.solve(b);
    return x;
}

Mat HODLR_Tree::solveNonLeafNonSPD(int j, int k, Mat b) 
{
    int r0 = tree[j][k]->rank[0];
    int r1 = tree[j][k]->rank[1];
    int n0 = tree[j][k]->c_size[0];
    int n1 = tree[j][k]->c_size[1];
    int r  = b.cols();

    // Initializing the temp matrix upon which inv(I + UV^T) is applied:
    //        |           |
    //        |  V1T * b  |
    //        |           |
    // temp = -------------
    //        |           |
    //        |  V0T * b  |
    //        |           |
    //        -------------

    Mat temp(r0 + r1, r);
    temp << tree[j][k]->V_factor[1].transpose() * b.block(n0, 0, n1, r),
            tree[j][k]->V_factor[0].transpose() * b.block(0,  0, n0, r);
    temp = tree[j][k]->K_factor_LU.solve(temp);
    
    Mat y(n0 + n1, r);
    
    // y = U * (1 + VT * U)^{-1} * VT * b
    y << tree[j][k]->U_factor[0] * temp.block(0,  0, r0, r), 
         tree[j][k]->U_factor[1] * temp.block(r0, 0, r1, r);
    
    // This is using Sherman Morrison Woodbury Formula:
    // x = b - U * (1 + VT * U)^{-1} * VT * b
    return(b - y);
}

Mat HODLR_Tree::solveNonSPD(Mat b) 
{   
    int start, size;
    Mat x = Mat::Zero(b.rows(),b.cols());
    
    int r = b.cols();

    // Solving over the leaf nodes:
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        start = tree[n_levels][k]->n_start;
        size  = tree[n_levels][k]->n_size;

        x.block(start, 0, size, r) = this->solveLeafNonSPD(k, b.block(start, 0, size, r));
    }

    b = x;
    
    // Solving over nonleaf levels:
    for(int j = n_levels - 1; j >= 0; j--) 
    {
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            start = tree[j][k]->n_start;
            size  = tree[j][k]->n_size;
            
            x.block(start, 0, size, r) = this->solveNonLeafNonSPD(j, k, b.block(start, 0, size, r));
        }

        b = x;
    }

    return x;
}

dtype HODLR_Tree::logDeterminantNonSPD()
{
    dtype log_det = 0.0;
    for(int j = n_levels; j >= 0; j--) 
    {
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            if(tree[j][k]->K.size() > 0)
            {
                for(int l = 0; l < tree[j][k]->K_factor_LU.matrixLU().rows(); l++) 
                {   
                    log_det += log(fabs(tree[j][k]->K_factor_LU.matrixLU()(l,l)));
                }
            }
        }
    }

    return(log_det);
}
