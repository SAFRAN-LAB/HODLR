#ifndef __HODLR_Node__
#define __HODLR_Node__

#include <Eigen/Dense>
#include "HODLR_Matrix.hpp"
#include "LowRank.hpp"

class HODLR_Node 
{
friend class HODLR;

// All methods are declared as private since all usage happens from
// the friend class HODLR_Tree:
private:
    HODLR_Node(int node_number, int level_number, int local_number, 
               int n_start, int n_size, double tolerance);
    
    // Storing the information passed to constructor as attribute:
    int node_number, level_number, local_number;
    int n_start, n_size;
    // Storing the start locations and sizes for the children of the node:
    int c_start[2], c_size[2];
    // Tolerance for the computation carried out:
    double tolerance;

    // This stores the matrix directly(i.e at the leaf level)
    Mat K;

    // Variables and methods needed for HODLR solver
    Mat U[2], V[2];
    Mat U_factor[2], V_factor[2];
    Eigen::PartialPivLU<Mat> K_factor_LU;
    int rank[2];
    // Variables needed for symmetric factorization
    // Unitary matrices:
    Mat Q[2];
    Eigen::LLT<Mat> K_factor_LLT;

    // Methods for Leaf Nodes:
    void assembleLeafNode(HODLR_Matrix* A);
    void matmatProductLeaf(Mat x, Mat& b);

    // Methods for Non-leaf Nodes:
    void assembleNonLeafNode(LowRank* A, bool is_sym);
    void matmatProductNonLeaf(Mat x, Mat& b);

    // Method to print the parameters of the node(mainly used to debug)
    void printNodeDetails();
};

#endif /*__HODLR_Node__*/
