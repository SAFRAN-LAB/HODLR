//
//  HODLR_Tree.hpp
//  
//
//  Created by Sivaram Ambikasaran on 11/8/13.
//
//

#ifndef __HODLR_TREE_HPP__
#define __HODLR_TREE_HPP__

#include <Eigen/Dense>
#include "HODLR_Node.hpp"

using namespace Eigen;

class HODLR_Tree {
private:
	friend class HODLR_Node;
	HODLR_Node* root;
	unsigned N;
	unsigned nLeaf;
	double lowRankTolerance;
	double determinant;
	double diagonal;

	void create_Tree(HODLR_Node*& node);
	void assemble_Matrix(HODLR_Node*& node);
	void matMatProduct(HODLR_Node*& node, MatrixXd& x, MatrixXd& b);
	void compute_Factor(HODLR_Node*& node);
	void set_Matrices_For_Inversion(HODLR_Node*& node);
	void solve(HODLR_Node*& node, MatrixXd& x);
	void compute_Determinant(HODLR_Node*& node);

public:
	HODLR_Tree(unsigned N, unsigned nLeaf);
	void assemble_Matrix(double diagonal, double lowRankTolerance);
	void matMatProduct(MatrixXd& x, MatrixXd& b);
	void compute_Factor();
	void solve(MatrixXd& b, MatrixXd& x);
	void compute_Determinant(double& determinant);
};

#endif /* defined(__HODLR_TREE_HPP__) */