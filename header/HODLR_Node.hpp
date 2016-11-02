#ifndef __HODLR_Node__
#define __HODLR_Node__

#include <Eigen/Dense>
#include <Eigen/QR>
#include <Eigen/Cholesky>
#include "HODLR_Matrix.hpp"

class HODLR_Node {
	friend class HODLR_Tree;
private:
	HODLR_Node(int nodeNumber, int levelNumber, int localNumber, int nStart, int nSize, double tolerance);
	int nodeNumber, levelNumber, localNumber;
	int nStart, nSize;
	int cStart[2], cSize[2];
	Eigen::MatrixXd U[2], V[2];
	Eigen::MatrixXd Ufactor[2], Vfactor[2];
	Eigen::MatrixXd K;
	Eigen::PartialPivLU<Eigen::MatrixXd> Kfactor;
	int rank[2];
	bool isLeaf;
	double tolerance;
	void assemble_Non_Leaf_Node(HODLR_Matrix* A);
	void assemble_Leaf_Node(HODLR_Matrix* A);
	void matmat_Product_Non_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b);
	void matmat_Product_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b);
	Eigen::MatrixXd Q[2];
	Eigen::MatrixXd Qfactor[2];
	Eigen::MatrixXd R;
	Eigen::LLT<Eigen::MatrixXd> llt;
	int sym_rank;
	void assemble_Symmetric_Non_Leaf_Node(HODLR_Matrix* A);
};

/********************************************************/
/*	PURPOSE OF EXISTENCE:	Constructor for the class	*/
/********************************************************/

/************/
/*	INPUTS	*/
/************/

/// nodeNumber      -   Node number of the HODLR Node in the tree.
/// levelNumber     -   Level number of the HODLR Node in the tree.
/// localNumber     -   Check for left or right child of the parent HODLR Node. '0' if the left child else '1'
/// nStart          -   Starting point of the HODLR node
/// nSize           -   Number of rows/columns of the HODLR node
/// tolerance       -   Permissible error in low-rank approximation of the off-diagonal blocks of the node

HODLR_Node::HODLR_Node(int nodeNumber, int levelNumber, int localNumber, int nStart, int nSize, double tolerance) {
	this->nodeNumber	=	nodeNumber;
	this->levelNumber	=	levelNumber;
	this->localNumber	=	localNumber;
	this->nStart		=	nStart;
	this->nSize			=	nSize;
	this->cStart[0]		=	nStart;
	this->cSize[0]		=	0.5*nSize;
	this->cStart[1]		=	nStart+cSize[0];
	this->cSize[1]		=	nSize-this->cSize[0];
	this->isLeaf		=	false;
	this->tolerance		=	tolerance;
}

/**********************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine to obtain and store low-rank approximation of the off-diagonal blocks of the node.*/
/**********************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// A  -  HODLR Matrix

void HODLR_Node::assemble_Non_Leaf_Node(HODLR_Matrix* A) {
	A->rook_Piv(cStart[0],cStart[1],cSize[0],cSize[1], tolerance, U[0], V[1], rank[0]);
	A->rook_Piv(cStart[1],cStart[0],cSize[1],cSize[0], tolerance, U[1], V[0], rank[1]);
}

/**************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine to store leaf node matrix blocks in K.*/
/**************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// A  -  HODLR Matrix

void HODLR_Node::assemble_Leaf_Node(HODLR_Matrix* A) {
	K	=	A->get_Matrix(nStart, nStart, nSize, nSize);
}

/**********************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	                                                                                          */
/**********************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/


void HODLR_Node::matmat_Product_Non_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	b.block(cStart[0],0,cSize[0],x.cols())+=(U[0]*(V[1].transpose()*x.block(cStart[1],0,cSize[1],x.cols())));
	b.block(cStart[1],0,cSize[1],x.cols())+=(U[1]*(V[0].transpose()*x.block(cStart[0],0,cSize[0],x.cols())));
}

/**********************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	                                                                                          */
/**********************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/


void HODLR_Node::matmat_Product_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	b.block(nStart,0,nSize,x.cols())+=K*x.block(nStart,0,nSize,x.cols());
}

/*************************************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine to obtain and store low-rank approximation of the off-diagonal block of the node for symmetric factorization.*/
/*************************************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// A  -  HODLR Matrix

void HODLR_Node::assemble_Symmetric_Non_Leaf_Node(HODLR_Matrix* A) {
	A->rook_Piv(cStart[0],cStart[1],cSize[0],cSize[1], tolerance, Q[0], Q[1], sym_rank);
}
#endif /*__HODLR_Node__*/
