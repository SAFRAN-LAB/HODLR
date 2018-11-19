#include "HODLR_Node.hpp"

HODLR_Node::HODLR_Node(int node_number, int level_number, int local_number, 
                       int n_start, int n_size, double tolerance
                      ) 
{
    // Storing the passed parameters as the attributes of the created object:
    this->node_number  = node_number;
    this->level_number = level_number;
    this->local_number = local_number;
    this->n_start      = n_start;
    this->n_size       = n_size;
    // Start location and size for the child on the left:
    this->c_start[0]   = n_start;
    this->c_size[0]    = 0.5 * n_size;
    // Start location and size for the child on the right:
    this->c_start[1]   = n_start + c_size[0];
    this->c_size[1]    = n_size  - c_size[0];
    // By default we set that created node isn't a leaf node:
    this->is_leaf      = false;
    this->tolerance    = tolerance;
}

void HODLR_Node::assembleLeafNode(HODLR_Matrix* A) 
{
    // At the leaf level we are just going to be building the matrix
    // directly since it's a full rank block:
    K = A->getMatrix(n_start, n_start, n_size, n_size);
}

void HODLR_Node::matmatProductLeaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) 
{
    b.block(n_start, 0, n_size, x.cols()) += K * x.block(n_start, 0, n_size, x.cols());
}

void HODLR_Node::assembleNonLeafNode(HODLR_Matrix* A) 
{
    A->rookPiv(c_start[0], c_start[1], c_size[0], c_size[1], tolerance, U[0], V[1], rank[0]);
    A->rookPiv(c_start[1], c_start[0], c_size[1], c_size[0], tolerance, U[1], V[0], rank[1]);
}

void HODLR_Node::matmatProductNonLeaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) 
{
    b.block(c_start[0], 0, c_size[0], x.cols()) += 
    (U[0] * (V[1].transpose() * x.block(c_start[1], 0, c_size[1], x.cols())));

    b.block(c_start[1], 0, c_size[1], x.cols()) += 
    (U[1] * (V[0].transpose() * x.block(c_start[0], 0, c_size[0], x.cols())));
}
