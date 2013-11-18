//
//  HODLR_Node.hpp
//  
//
//  Created by Sivaram Ambikasaran on 11/8/13.
//
//

#ifndef __HODLR_NODE_HPP__
#define __HODLR_NODE_HPP__

#include <Eigen/Dense>

using namespace Eigen;

class HODLR_Node {
	friend class HODLR_Tree;
private:
	HODLR_Node* parent;
	HODLR_Node* child[2];

//	Variables of interest for both leaf and non-leaf;
	unsigned levelNumber;	//	Level number of the node;
	unsigned nodeNumber;	//	Node number is either 0 or 1;
	unsigned nStart;	//	nStart is the starting index of node;
	unsigned nSize;		//	nSize is the size of the node;
	MatrixXd K;		//	At leaf level, stores the self interaction; At non-leaf stores the matrix [I, V1inverse*U1inverse; V0inverse*U0inverse, I];
	FullPivLU<MatrixXd> Kinverse;	//	Stores Factorization of K;
	double determinant;	//	Stores K.determinant();

//	Variables of interest for non-leaf;
	unsigned nRank[2];	//	nRank[0] rank of K01; nRank[1] rank of K10;
	MatrixXd U[2];		//	Column basis of low-rank interaction of children at non-leaf;
	MatrixXd V[2];		//	Row basis of low-rank interaction ofchildren at non-leaf;
	MatrixXd Uinverse[2];	//	Column basis of low-rank interaction of children of inverse at non-leaf.
	MatrixXd Vinverse[2];	//	Row basis of low-rank interaction of children of inverse at non-leaf.

//	Variables of interest at leaf;
	bool isLeaf;		//	If node is a leaf, it takes the value TRUE;

	HODLR_Node(unsigned levelNumber, unsigned nodeNumber, unsigned nStart, unsigned nSize);

	void assemble_Matrices(double lowRankTolerance, VectorXd& diagonal);

	void matrix_Matrix_Product(MatrixXd& x, MatrixXd& b);

	void set_UV_Inversion();

	void compute_K();

	void compute_Inverse();

	void apply_Inverse(MatrixXd& matrix, unsigned mStart);

	void compute_Determinant();

	~HODLR_Node();
};

#endif /* defined(__HODLR_NODE_HPP__) */