#ifndef __HODLR_Tree__
#define __HODLR_Tree__

#include <Eigen/Dense>
#include <vector>

#include "HODLR_Matrix.hpp"
#include "HODLR_Node.hpp"

class HODLR_Tree 
{
private:
    int N;
    int n_levels;
    double tolerance;
    std::vector<int> nodes_in_level;
    HODLR_Matrix* A;
    bool is_sym, is_pd;
    
    // Vector of levels(which contain nodes) thereby giving the tree:
    std::vector<std::vector<HODLR_Node*>> tree;
    void createTree();
    void createRoot();
    void createChildren(int level_number, int node_number);

    // Methods needed for the nonSPD HODLR solver:
    void factorizeLeafNonSPD(int node_number);
    void factorizeNonLeafNonSPD(int level_number, int node_number);
    void factorizeNonSPD();
    MatrixXd solveLeafNonSPD(int node_number, MatrixXd b);
    MatrixXd solveNonLeafNonSPD(int level_number, int node_number, MatrixXd b);
    MatrixXd solveNonSPD(MatrixXd b);
    double logDeterminantNonSPD();

    // Methods needed for the SPD HODLR solver:
    void factorizeLeafSPD(int node_number);
    void factorizeNonLeafSPD(int level_number, int node_number);
    void factorizeSPD();
    void qr(int level_number, int node_number);
    void qrForLevel(int level_number);
    MatrixXd solveLeafSymmetricFactor(int node_number, MatrixXd b);
    MatrixXd solveNonLeafSymmetricFactor(int level_number, int node_number, MatrixXd b);
    MatrixXd solveSymmetricFactor(MatrixXd b);
    MatrixXd solveLeafSymmetricFactorTranspose(int node_number, MatrixXd b);
    MatrixXd solveNonLeafSymmetricFactorTranspose(int level_number, int node_number, MatrixXd b);
    MatrixXd solveSymmetricFactorTranspose(MatrixXd b);
    MatrixXd solveSPD(MatrixXd b);
    MatrixXd SymmetricFactorNonLeafProduct(int level_number, int node_number, MatrixXd b);
    MatrixXd SymmetricFactorTransposeNonLeafProduct(int level_number, int node_number, MatrixXd b);

    double logDeterminantSPD();

public:
    HODLR_Tree(int n_levels, double tolerance, HODLR_Matrix* A);
    ~HODLR_Tree();

    //  Methods for HODLR solver
    void assembleTree(bool is_sym = false, bool is_pd = false);
    // Gives the box details of the prescribed box and level number:
    void printNodeDetails(int level_number, int node_number);
    // Lists details of all boxes in the tree
    void printTreeDetails();
    void plotTree();
    MatrixXd matmatProduct(MatrixXd x);
    void factorize();
    MatrixXd solve(MatrixXd b);
    MatrixXd symmetricFactorProduct(MatrixXd x);
    MatrixXd symmetricFactorTransposeProduct(MatrixXd x);
    MatrixXd getSymmetricFactor();
    double logDeterminant();
};

#endif /*__HODLR_Tree__*/
