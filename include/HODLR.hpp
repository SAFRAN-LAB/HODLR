#ifndef __HODLR__
#define __HODLR__

#include <Eigen/Dense>
#include "LowRank.hpp"
#include "HODLR_Node.hpp"

class HODLR 
{
private:
    int n_levels;
    double tolerance;
    LowRank* F;
    std::vector<int> nodes_in_level;
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
    Mat solveLeafNonSPD(int node_number, Mat b);
    Mat solveNonLeafNonSPD(int level_number, int node_number, Mat b);
    Mat solveNonSPD(Mat b);
    dtype logDeterminantNonSPD();

    // Methods needed for the SPD HODLR solver:
    void factorizeLeafSPD(int node_number);
    void factorizeNonLeafSPD(int level_number, int node_number);
    void factorizeSPD();
    void qr(int level_number, int node_number);
    void qrForLevel(int level_number);
    Mat solveLeafSymmetricFactor(int node_number, Mat b);
    Mat solveNonLeafSymmetricFactor(int level_number, int node_number, Mat b);
    Mat solveSymmetricFactor(Mat b);
    Mat solveLeafSymmetricFactorTranspose(int node_number, Mat b);
    Mat solveNonLeafSymmetricFactorTranspose(int level_number, int node_number, Mat b);
    Mat solveSymmetricFactorTranspose(Mat b);
    Mat solveSPD(Mat b);
    Mat SymmetricFactorNonLeafProduct(int level_number, int node_number, Mat b);
    Mat SymmetricFactorTransposeNonLeafProduct(int level_number, int node_number, Mat b);
    dtype logDeterminantSPD();

public:
    
    // Size of the matrix considered:
    int N;

    HODLR(int N, int M, double tolerance);
    
    ~HODLR();

    //  Methods for HODLR solver
    // Gives the box details of the prescribed box and level number:
    void printNodeDetails(int level_number, int node_number);
    // Lists details of all boxes in the tree
    void printTreeDetails();
    void plotTree(std::string image_name);
    void assemble(HODLR_Matrix* A, std::string lowrank_type = "rookPivoting", 
                  bool is_sym = false, bool is_pd = false
                 );
    Mat matmatProduct(Mat x);
    void factorize();
    Mat solve(Mat b);
    Mat symmetricFactorProduct(Mat x);
    Mat symmetricFactorTransposeProduct(Mat x);
    Mat getSymmetricFactor();
    dtype logDeterminant();
};

#endif /*__HODLR__*/
