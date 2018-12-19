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
    bool is_spd, is_sym;
	
    // Vector of levels(which contain nodes) thereby giving the tree:
    std::vector<std::vector<HODLR_Node*>> tree;
	void createTree();
	void createRoot();
	void createChildren(int level_number, int node_number);

	// Variables and methods needed for HODLR solver
	void factorizeLeaf(int node_number);
	void factorizeNonLeaf(int level_number, int node_number);
    void qr(int level_number, int node_number);
    void qrForLevel(int level_number);
	MatrixXd solveLeaf(int node_number, MatrixXd b);
	MatrixXd solveNonLeaf(int level_number, int node_number, MatrixXd b);

public:
	HODLR_Tree(int n_levels, double tolerance, HODLR_Matrix* A);
	~HODLR_Tree();

	//  Methods for HODLR solver
	void assembleTree(bool is_spd = false, bool is_sym = false);
    // Gives the box details of the prescribed box and level number:
    void printNodeDetails(int level_number, int box_number);
    // Lists details of all boxes in the tree
    void printTreeDetails();
	void factorize();
	void matmatProduct(MatrixXd x, MatrixXd& b);
    void symmetricFactorProduct(MatrixXd x, MatrixXd& b);
    void symmetricFactorTransposeProduct(MatrixXd x, MatrixXd& b);
	double logDeterminant();
	MatrixXd solve(MatrixXd b);
    void plotTree();
};

#endif /*__HODLR_Tree__*/
