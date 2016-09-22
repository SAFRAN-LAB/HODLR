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
	//Symmetric factorization nodes
	Eigen::MatrixXd Q[2];
	Eigen::MatrixXd Qfactor[2];
	Eigen::MatrixXd R;
	Eigen::MatrixXd Rfactor[2];
	Eigen::MatrixXd X;
	Eigen::LLT<Eigen::MatrixXd> llt;
	int sym_rank;
	void assemble_Symmetric_Non_Leaf_Node(HODLR_Matrix* A);
};

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

void HODLR_Node::assemble_Non_Leaf_Node(HODLR_Matrix* A) {
	// std::cout << "\nStart assemble_Non_Leaf_Node\n";
	A->rook_Piv(cStart[0],cStart[1],cSize[0],cSize[1], tolerance, U[0], V[1], rank[0]);
	A->rook_Piv(cStart[1],cStart[0],cSize[1],cSize[0], tolerance, U[1], V[0], rank[1]);
	// std::cout << "\nDone assemble_Non_Leaf_Node\n";
}

void HODLR_Node::assemble_Leaf_Node(HODLR_Matrix* A) {
	// std::cout << "\nStart assemble_Leaf_Node\n";
	K	=	A->get_Matrix(nStart, nStart, nSize, nSize);
	// std::cout << "\nDone assemble_Leaf_Node\n";
}

void HODLR_Node::matmat_Product_Non_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	// std::cout << "\nStart matmat_Product_Non_Leaf\n";
	// std::cout << U[0].rows() << "\n";
	// std::cout << U[0].cols() << "\n";
	// std::cout << V[1].rows() << "\n";
	// std::cout << V[1].cols() << "\n";
	// std::cout << cSize[1] << "\n";
	b.block(cStart[0],0,cSize[0],x.cols())+=(U[0]*(V[1].transpose()*x.block(cStart[1],0,cSize[1],x.cols())));
	b.block(cStart[1],0,cSize[1],x.cols())+=(U[1]*(V[0].transpose()*x.block(cStart[0],0,cSize[0],x.cols())));
	// std::cout << "\nDone matmat_Product_Non_Leaf\n";
}

void HODLR_Node::matmat_Product_Leaf(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	// std::cout << "\nStart matmat_Product_Leaf\n";
	b.block(nStart,0,nSize,x.cols())+=K*x.block(nStart,0,nSize,x.cols());
	// std::cout << "\nDone matmat_Product_Leaf\n";
}

void HODLR_Node::assemble_Symmetric_Non_Leaf_Node(HODLR_Matrix* A) {
	// std::cout << "\nStart assemble_Non_Leaf_Node\n";
	/*Eigen::MatrixXd USym[2];
	Eigen::MatrixXd R[2];
	A->rook_Piv(cStart[0],cStart[1],cSize[0],cSize[1], tolerance, USym[0], USym[1], sym_rank);

	int minmn[2];
	minmn[0] = std::min(USym[0].rows(), USym[0].cols());
	minmn[1] = std::min(USym[1].rows(), USym[1].cols());
	
	HouseholderQR<MatrixXd> qr(USym[0]);
	Eigen::MatrixXd::Identity(USym[0].rows(), minmn[0]);
    	A->Q[0] = qr.householderQ()*(Eigen::MatrixXd::Identity((USym[0].rows(), minmn[0]));
	R[0] = qr.matrixQR().block(0,0,minmn[0],USym[0].cols()).triangularView<Eigen::Upper>();
    	
	qr(USym[1]);
  	Eigen::MatrixXd::Identity(USym[1].rows(), minmn[1]);
    	A->Q[1] = qr.householderQ()*(Eigen::MatrixXd::Identity((USym[1].rows(), minmn[1]));
	R[1] = qr.matrixQR().block(0,0,minmn[1],USym[1].cols()).triangularView<Eigen::Upper>();
  	
	A->R = R[0]*R[1]*transpose();*/

	A->rook_Piv(cStart[0],cStart[1],cSize[0],cSize[1], tolerance, Q[0], Q[1], sym_rank);


	// std::cout << "\nDone assemble_Non_Leaf_Node\n";
}
#endif /*__HODLR_Node__*/
