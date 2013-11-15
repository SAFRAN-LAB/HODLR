//
//  HODLR_Node.cpp
//  
//
//  Created by Sivaram Ambikasaran on 11/8/13.
//
//

#include <Eigen/Dense>
#include <iostream>

#include "HODLR_Node.hpp"
#include "get_Matrix.hpp"
#include "partial_Piv_LU.hpp"

using namespace std;
using namespace Eigen;

HODLR_Node::HODLR_Node(unsigned levelNumber, unsigned nodeNumber, unsigned nStart, unsigned nSize) {
	this->levelNumber	=	levelNumber;
	this->nodeNumber	=	nodeNumber;
	this->nStart			=	nStart;
	this->nSize				=	nSize;
	this->parent			=	NULL;
	this->child[0]		=	NULL;
	this->child[1]		=	NULL;
}

void HODLR_Node::assemble_Matrices(double lowRankTolerance, double diagonal) {
	if (isLeaf	==	true) {
		get_Matrix(nStart, nStart, nSize, nSize, K);
		for (unsigned k=0; k<nSize; ++k) {
			K(k,k)	=	diagonal;
		}
	}
	else if (isLeaf	==	false) {

		partial_Piv_LU(child[0]->nStart, child[1]->nStart, child[0]->nSize, child[1]->nSize, lowRankTolerance, nRank[0], U[0], V[1]);

		partial_Piv_LU(child[1]->nStart, child[0]->nStart, child[1]->nSize, child[0]->nSize, lowRankTolerance, nRank[1], U[1], V[0]);

	}
}

void HODLR_Node::matrix_Matrix_Product(MatrixXd& x, MatrixXd& b) {
	unsigned n	=	x.cols();

	if (isLeaf	==	true) {
		b.block(nStart, 0, nSize, n)	=	b.block(nStart, 0, nSize, n)	+	K*x.block(nStart, 0, nSize, n);
	}
	else if (isLeaf	==	false) {
		unsigned n0	=	child[0]->nStart;
		unsigned n1	=	child[1]->nStart;
		unsigned m0	=	child[0]->nSize;
		unsigned m1	=	child[1]->nSize;
		b.block(n0, 0, m0, n)	=	b.block(n0, 0, m0, n)	+	U[0]*(V[1]*x.block(n1, 0, m1, n));
		b.block(n1, 0, m1, n)	=	b.block(n1, 0, m1, n)	+	U[1]*(V[0]*x.block(n0, 0, m0, n));
	}
}

void HODLR_Node::set_UV_Inversion() {
	for (unsigned k=0; k<2; ++k) {
		Uinverse[k]	=	U[k];
		Vinverse[k]	=	V[k];
	}
}

void HODLR_Node::compute_K() {
	if (isLeaf	==	false) {
		unsigned m0	=	V[0].rows();
		unsigned m1	=	V[1].rows();
		K	=	MatrixXd::Identity(m0+m1, m0+m1);

		K.block(0, m1, m1, m0)	=	Vinverse[1]*Uinverse[1];
		K.block(m1, 0, m0, m1)	=	Vinverse[0]*Uinverse[0];
	}
}

void HODLR_Node::compute_Inverse() {
	Kinverse.compute(K);
}

void HODLR_Node::apply_Inverse(MatrixXd& matrix, unsigned mStart) {
	unsigned n	=	matrix.cols();
	unsigned start	=	nStart-mStart;
	if (isLeaf	==	true) {
		matrix.block(start, 0, nSize, n)	=	Kinverse.solve(matrix.block(start, 0, nSize, n));
	}
	else if (isLeaf	==	false) {
		//	Computes temp		=	Vinverse*matrix

		MatrixXd temp(nRank[0]+nRank[1], n);

		temp.block(0, 0, nRank[0] , n)			=	Vinverse[1]*matrix.block(start+child[0]->nSize, 0 , child[1]->nSize, n);

		temp.block(nRank[0], 0, nRank[1] , n)	=	Vinverse[0]*matrix.block(start, 0 , child[0]->nSize, n);
		
		//	Computes tempSolve	=	Kinverse\temp

		MatrixXd tempSolve	=	Kinverse.solve(temp);
		
		//	Computes matrix		=	matrix-Uinverse*tempSolve

		matrix.block(start, 0, child[0]->nSize, n)					=	matrix.block(start, 0, child[0]->nSize, n)	-	Uinverse[0]*tempSolve.block(0, 0, nRank[0], n);
		matrix.block(start + child[0]->nSize, 0, child[1]->nSize, n)=	matrix.block(start + child[0]->nSize, 0, child[1]->nSize, n)	-	Uinverse[1]*tempSolve.block(nRank[0], 0, nRank[1], n);
	}
}

void HODLR_Node::compute_Determinant() {
	determinant	=	log(fabs(K.determinant()));
}

HODLR_Node::~HODLR_Node() {
	delete child[0];
	delete child[1];
	child[0]	=	NULL;
	child[1]	=	NULL;
}
