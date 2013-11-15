//
//  HODLR_Tree.cpp
//  
//
//  Created by Sivaram Ambikasaran on 11/8/13.
//
//

#include <Eigen/Dense>
#include "HODLR_Tree.hpp"
#include "get_Matrix.hpp"
#include <iostream>

using namespace Eigen;

HODLR_Tree::HODLR_Tree(unsigned N, unsigned nLeaf) {
	this->N		=	N;
	this->nLeaf	=	nLeaf;
	root		=	new HODLR_Node(0, 0, 0, N);
	create_Tree(root);
}

void HODLR_Tree::create_Tree(HODLR_Node*& node) {
	if (node->nSize <= nLeaf) {
		node->isLeaf	=	true;
	}
	else {
		node->isLeaf	=	false;
		unsigned n0		=	node->nSize/2;
		unsigned n1		=	node->nSize-n0;
		node->child[0]	=	new HODLR_Node(node->levelNumber+1, 0, node->nStart, n0);
		node->child[0]->parent	=	node;
		node->child[1]	=	new HODLR_Node(node->levelNumber+1, 1, node->nStart+n0, n1);
		node->child[1]->parent	=	node;
		create_Tree(node->child[0]);
		create_Tree(node->child[1]);
	}
}

void HODLR_Tree::assemble_Matrix(double diagonal, double lowRankTolerance) {
	this->lowRankTolerance	=	lowRankTolerance;
	this->diagonal		=	diagonal;
	set_Locations(N);
	assemble_Matrix(root);
}

void HODLR_Tree::assemble_Matrix(HODLR_Node*& node) {
	if (node) {
		node->assemble_Matrices(lowRankTolerance, diagonal);
		assemble_Matrix(node->child[0]);
		assemble_Matrix(node->child[1]);
	}
}

void HODLR_Tree::matMatProduct(MatrixXd& x, MatrixXd& b) {

	b	=	MatrixXd::Zero(N, x.cols());
	matMatProduct(root, x, b);
}

void HODLR_Tree::matMatProduct(HODLR_Node*& node, MatrixXd& x, MatrixXd& b) {

	if (node) {
		node->matrix_Matrix_Product(x, b);
		matMatProduct(node->child[0], x, b);
		matMatProduct(node->child[1], x, b);
	}
}

void HODLR_Tree::compute_Factor() {
	set_Matrices_For_Inversion(root);
	compute_Factor(root);
}

void HODLR_Tree::compute_Factor(HODLR_Node*& node) {
	if (node) {
		for (unsigned k=0; k<2; ++k) {
			compute_Factor(node->child[k]);
		}
		node->compute_K();
		node->compute_Inverse();
		HODLR_Node*& mynode	=	node->parent;
		unsigned number		=	node->nodeNumber;
		unsigned mStart		=	node->nStart;
		while (mynode) {
			node->apply_Inverse(mynode->Uinverse[number], mStart);
			number		=	mynode->nodeNumber;
			mStart		=	mynode->nStart;
			mynode		=	mynode->parent;
		}
	}
}

void HODLR_Tree::set_Matrices_For_Inversion(HODLR_Node*& node) {
	if (node) {
		node->set_UV_Inversion();
		set_Matrices_For_Inversion(node->child[0]);
		set_Matrices_For_Inversion(node->child[1]);
	}
}

void HODLR_Tree::solve(MatrixXd& b, MatrixXd& x) {
	x	=	b;
	solve(root, x);
	
}

void HODLR_Tree::solve(HODLR_Node*& node, MatrixXd& x) {
	if (node) {
		solve(node->child[0], x);
		solve(node->child[1], x);
		node->apply_Inverse(x, 0);
	}
}

void HODLR_Tree::compute_Determinant(double& determinant) {
	this->determinant	=	0;
	compute_Determinant(root);
	determinant		=	this->determinant;
}

void HODLR_Tree::compute_Determinant(HODLR_Node*& node) {
	if (node) {
		for (unsigned k=0; k<2; ++k) {
			compute_Determinant(node->child[k]);
			node->compute_Determinant();
		}
		determinant	=	determinant+node->determinant;
	}
}