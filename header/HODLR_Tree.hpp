#ifndef __HODLR_Tree__
#define __HODLR_Tree__

#include <Eigen/Dense>
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <vector>

#include "HODLR_Matrix.hpp"
#include "HODLR_Node.hpp"

class HODLR_Tree {
	private:
		int N;
		int nLevels;
		double tolerance;
		std::vector<int> nodesInLevel;
		HODLR_Matrix* A;
		std::vector< std::vector<HODLR_Node*> > tree;
		void createTree();
		void createRoot();
		void createChildren(int j, int k);


		//  Variables and methods needed for HODLR solver
		void factorize_Leaf(int k);
		void factorize_Non_Leaf(int j, int k);
		Eigen::MatrixXd solve_Leaf(int k, Eigen::MatrixXd b);
		Eigen::MatrixXd solve_Non_Leaf(int j, int k, Eigen::MatrixXd b);


		//  Variables and methods needed for symmetric factorization
		void qr(int j, int k);
		void qr_Level(int level);
		void factorize_Symmetric_Leaf(int k);
		void factorize_Symmetric_Non_Leaf(int j, int k);
		Eigen::MatrixXd solve_Symmetric_Factor_Leaf(int k, Eigen::MatrixXd b);
		Eigen::MatrixXd solve_Symmetric_Factor_Non_Leaf(int j, int k,  Eigen::MatrixXd b);

		void small_Cholesky(int j, int k);

		Eigen::MatrixXd mult_Symmetric_Factor_Non_Leaf(int j, int k, Eigen::MatrixXd b);
		Eigen::MatrixXd mult_Symmetric_Factor_transpose_Non_Leaf(int j, int k, Eigen::MatrixXd b);

		Eigen::MatrixXd solve_Symmetric_Factor_Transpose_Non_Leaf(int j, int k, Eigen::MatrixXd b);
		Eigen::MatrixXd solve_Symmetric_Factor_Transpose_Leaf(int k, Eigen::MatrixXd b);

	public:
		HODLR_Tree(int nLevels, double tolerance, HODLR_Matrix* A);
		~HODLR_Tree();
		virtual double get_Matrix_Element(int j, int k) {
			return 0.0;
		}

		//  Methods for HODLR solver
		void assemble_Tree();
		void factorize();
		double determinant();
		void matmat_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b);
		Eigen::MatrixXd solve(Eigen::MatrixXd b);

		//  Methods for symmetric factorization
		void assembleSymmetricTree();
		void symmetric_factorize();
		double symmetric_Determinant();
		void matmat_Symmetric_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b);
		Eigen::MatrixXd symmetric_Solve(Eigen::MatrixXd b);
		Eigen::MatrixXd symmetric_Factor_Solve(Eigen::MatrixXd b);
		Eigen::MatrixXd symmetric_Factor_Transpose_Solve(Eigen::MatrixXd b);
		Eigen::MatrixXd symmetric_Factor_Product(Eigen::MatrixXd b);
		Eigen::MatrixXd symmetric_Factor_Transpose_Product(Eigen::MatrixXd b);

		//  Don't use this for large N; Essentially builds the entire symmetric factor; Only for testing purposes;
		Eigen::MatrixXd build_Symmetric_Matrix_Factor();
};

/*******************************************************/
/*	PURPOSE OF EXISTENCE: Constructor for the class    */
/*******************************************************/

/************/
/*	INPUTS	*/
/************/

/// nLevels      -   Number of levels in the HODLR Tree
/// tolerance    -   Relative error for approximating the off-diagonal blocks
/// A            -   HODLR matrix for the kernel matrix

HODLR_Tree::HODLR_Tree(int nLevels, double tolerance, HODLR_Matrix* A) {
	// // std::cout << "\nStart HODLR_Tree\n";
	this->nLevels	=	nLevels;
	this->tolerance	=	tolerance;
	this->A			=	A;
	this->N			=	A->N;
	nodesInLevel.push_back(1);
	for (int j=1; j<=nLevels; ++j) {
		nodesInLevel.push_back(2*nodesInLevel.back());
	}
	// // std::cout << "\nDone HODLR_Tree\n";
	createTree();
}

/*******************************************************/
/*	PURPOSE OF EXISTENCE: Builds the root of HODLR Tree*/
/*******************************************************/

void HODLR_Tree::createRoot() {
	// // std::cout << "\nStart createRoot\n";
	HODLR_Node* root	=	new HODLR_Node(0, 0, 0, 0, N, tolerance);
	std::vector<HODLR_Node*> level;
	level.push_back(root);
	tree.push_back(level);
	// // std::cout << "\nDone createRoot\n";
}

/******************************************************************************/
/*	PURPOSE OF EXISTENCE: Builds HODLR child nodes for HODLR nodes in the tree*/
/******************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number

void HODLR_Tree::createChildren(int j, int k) {
	// // std::cout << "\nStart createChildren\n";
	//	Adding left child
	HODLR_Node* left	=	new HODLR_Node(2*j, k+1, 0, tree[j][k]->cStart[0], tree[j][k]->cSize[0], tolerance);
	tree[j+1].push_back(left);

	//	Adding right child
	HODLR_Node* right	=	new HODLR_Node(2*j+1, k+1, 1, tree[j][k]->cStart[1], tree[j][k]->cSize[1], tolerance);
	tree[j+1].push_back(right);
	// // std::cout << "\nDone createChildren\n";
}

/****************************************************************/
/*	PURPOSE OF EXISTENCE: Builds HODLR nodes for HODLR Matrix A */
/****************************************************************/

void HODLR_Tree::createTree() {
	// // std::cout << "\nStart createTree\n";
	createRoot();
	for (int j=0; j<nLevels; ++j) {
		std::vector<HODLR_Node*> level;
		tree.push_back(level);
		for (int k=0; k<nodesInLevel[j]; ++k) {
			createChildren(j,k);
		}
	}
	// // std::cout << "\nDone createTree\n";
}

/************************************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE: Obtains a factorization of the leaf nodes and computes the low rank approximations of the off-diagonal blocks, Z=UV'. */
/************************************************************************************************************************************************/

void HODLR_Tree::assemble_Tree() {
	// // std::cout << "\nStart assemble_Tree\n";
	for (int j=0; j<nLevels; ++j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->assemble_Non_Leaf_Node(A);
		}
	}
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->assemble_Leaf_Node(A);
	}
	// // std::cout << "\nDone assemble_Tree\n";
}

/*******************************************************/
/*	PURPOSE OF EXISTENCE:	 Matrix-matrix product     */
/*******************************************************/

/************/
/*	INPUTS	*/
/************/

/// x   -   Matrix to be multiplied on the right of the HODLR matrix

/************/
/*	OUTPUTS	*/
/************/

/// b   -   Matrix matrix product

void HODLR_Tree::matmat_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	// // std::cout << "\nStart matmat_Product\n";
	int r	=	x.cols();
	b		=	Eigen::MatrixXd::Zero(N,r);
	for (int j=0; j<nLevels; ++j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->matmat_Product_Non_Leaf(x, b);
		}
	}
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->matmat_Product_Leaf(x, b);
	}
	// // std::cout << "\nDone matmat_Product\n";
}

//	Factorization begins from here

/****************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining factorisations of leaf nodes	*/
/****************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// k           -   Leaf Node number

void HODLR_Tree::factorize_Leaf(int k) {
	// std::cout << "\nStart factorize Leaf: " << k << "\n";
	tree[nLevels][k]->Kfactor.compute(tree[nLevels][k]->K);
	int parent	=	k;
	int child	=	k;
	int size	=	tree[nLevels][k]->nSize;
	int tstart, r;
    #pragma omp parallel for
	for (int l=nLevels-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[nLevels][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->rank[child];
		tree[l][parent]->Ufactor[child].block(tstart,0,size,r)	=	solve_Leaf(k, tree[l][parent]->Ufactor[child].block(tstart,0,size,r));
	}
}

/************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving Kx = b, where K is the leaf node k. */
/************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of K multiplied with b


Eigen::MatrixXd HODLR_Tree::solve_Leaf(int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart solve Leaf: " << k << "\n";
	Eigen::MatrixXd x	=	tree[nLevels][k]->Kfactor.solve(b);
	// std::cout << "\nDone solve Leaf: " << k << "\n";
	return x;
}

/********************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining factorisations of non-leaf nodes	*/
/********************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number

void HODLR_Tree::factorize_Non_Leaf(int j, int k) {
	// std::cout << "\nStart factorize; Level: " << j << "Node: " << k << "\n";
	int r0	=	tree[j][k]->rank[0];
	int r1	=	tree[j][k]->rank[1];
	tree[j][k]->K	=	Eigen::MatrixXd::Identity(r0+r1, r0+r1);
	// std::cout << tree[j][k]->K.rows() << "\t" << tree[j][k]->K.cols() << "\n";
	tree[j][k]->K.block(0, r0, r0, r1)	=	tree[j][k]->Vfactor[1].transpose()*tree[j][k]->Ufactor[1];
	tree[j][k]->K.block(r0, 0, r1, r0)	=	tree[j][k]->Vfactor[0].transpose()*tree[j][k]->Ufactor[0];
	tree[j][k]->Kfactor.compute(tree[j][k]->K);
	int parent	=	k;
	int child	=	k;
	int size	=	tree[j][k]->nSize;
	int tstart, r;

	for (int l=j-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[j][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->rank[child];
		// std::cout << "\n" << size << "\n";
		tree[l][parent]->Ufactor[child].block(tstart,0,size,r)	=	solve_Non_Leaf(j, k, tree[l][parent]->Ufactor[child].block(tstart,0,size,r));
	}
	// std::cout << "\nDone factorize; Level: " << j << "Node: " << k << "\n";
}

/********************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving (I+UKV')x = b. The method uses Sherman-Morrison-Woodsbury formula to obtain x.	*/
/********************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	(I-U(inverse of (I+K))KV') multiplied by b

Eigen::MatrixXd HODLR_Tree::solve_Non_Leaf(int j, int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart Solve non leaf; Level: " << j << "Node: " << k << "\n";
	int r0	=	tree[j][k]->rank[0];
	int r1	=	tree[j][k]->rank[1];
	int n0	=	tree[j][k]->cSize[0];
	int n1	=	tree[j][k]->cSize[1];
	int r	=	b.cols();
	Eigen::MatrixXd temp(r0+r1, r);
	temp << tree[j][k]->Vfactor[1].transpose()*b.block(n0,0,n1,r),
	     tree[j][k]->Vfactor[0].transpose()*b.block(0,0,n0,r);
	temp	=	tree[j][k]->Kfactor.solve(temp);
	Eigen::MatrixXd y(n0+n1, r);
	y << tree[j][k]->Ufactor[0]*temp.block(0,0,r0,r), tree[j][k]->Ufactor[1]*temp.block(r0,0,r1,r);
	return (b-y);
}

/*********************************************************/
/*	PURPOSE OF EXISTENCE: Factorises the kernel matrix A.*/
/*********************************************************/

void HODLR_Tree::factorize() {
	// std::cout << "\nStart factorize...\n";
	//	Set things for factorization
	// std::cout << "Number of levels: " << nLevels << "\n";
	for (int j=0; j<=nLevels; ++j) {
    #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			for (int l=0; l<2; ++l) {
				tree[j][k]->Ufactor[l]	=	tree[j][k]->U[l];
				tree[j][k]->Vfactor[l]	=	tree[j][k]->V[l];
				// std::cout << "Level " << j << "; Node " << k << "; Matrix Size" << tree[j][k]->Ufactor[l].rows() << "\t" << tree[j][k]->Ufactor[l].cols() << "\n";
			}
		}
	}

	//	Factorize the leaf nodes
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		factorize_Leaf(k);
	}
	//	Factorize the non-leaf nodes
	for (int j=nLevels-1; j>=0; --j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			factorize_Non_Leaf(j, k);
		}
	}
	// std::cout << "\nEnd factorize...\n";
}

/**********************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Solves for the 	linear system Ax=b, where A is the kernel matrix. */
/**********************************************************************************************/

/************/
/*	INPUTS	*/
/************/

///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of kernel matrix multiplied by input matrix

Eigen::MatrixXd HODLR_Tree::solve(Eigen::MatrixXd b) {
	// std::cout << "\nStart solve...\n";
	int start, size;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();
	// #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r)	=	solve_Leaf(k, b.block(start, 0, size, r));
	}
	b=x;
	for (int j=nLevels-1; j>=0; --j) {
		// #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	solve_Non_Leaf(j, k, b.block(start, 0, size, r));
		}
		b=x;
	}
	// std::cout << "\nEnd solve...\n";
	return x;
}

/**********************************************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE: Obtains a symmetric factorization of the leaf nodes and computes the low rank approximations of the off-diagonal blocks, Z=UV' .*/
/**********************************************************************************************************************************************************/

void HODLR_Tree::assembleSymmetricTree() {
	std::cout << "\nStart Symmetric assemble_Tree\n";
	for (int j=0; j<nLevels; ++j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->assemble_Symmetric_Non_Leaf_Node(A);
		}
	}
	#pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->assemble_Leaf_Node(A);
		//Storing Cholesky decomposition of leaf nodes
		tree[nLevels][k]->llt.compute(tree[nLevels][k]->K);

	}
	std::cout << "\nDone Symmetric assemble_Tree\n";
}

/*******************************************************/
/*	PURPOSE OF EXISTENCE:	 Matrix-matrix product     */
/*******************************************************/

/************/
/*	INPUTS	*/
/************/

/// x   -   Matrix to be multiplied on the right of the HODLR matrix

/************/
/*	OUTPUTS	*/
/************/

/// b   -   Matrix matrix product

void HODLR_Tree::matmat_Symmetric_Product(Eigen::MatrixXd x, Eigen::MatrixXd& b) {
	// // std::cout << "\nStart matmat_Product\n";
	int r	=	x.cols();
	b		=	Eigen::MatrixXd::Zero(N,r);
	for (int j=0; j<nLevels; ++j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->matmat_Symmetric_Product_Non_Leaf(x, b);
		}
	}
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		tree[nLevels][k]->matmat_Product_Leaf(x, b);
	}
	// // std::cout << "\nDone matmat_Product\n";
}

/******************************************************************************************************************/
/*	PURPOSE OF EXISTENCE: Does a symmetric factorization of the symmetric positive-definite kernel matrix A = WW'.*/
/******************************************************************************************************************/

void HODLR_Tree::symmetric_factorize() {

	std::cout << "\nStart Symmetric factorize\n";

	for (int j=0; j<nLevels; ++j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			tree[j][k]->Qfactor[0]	=	tree[j][k]->Q[0];
			tree[j][k]->Qfactor[1]	=	tree[j][k]->Q[1];
			tree[j][k]->K = Eigen::MatrixXd::Identity(tree[j][k]->sym_rank, tree[j][k]->sym_rank);
		}
	}

	//	Factorize the leaf nodes
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		factorize_Symmetric_Leaf(k);
	}
	//	std::cout << "\nStart Leaf qr...\n";
	qr_Level(nLevels-1);
	//	std::cout << "\nEnd Leaf qr...\n";

	//	Factorize the non-leaf nodes
	for (int j=nLevels-1; j>0; --j) {
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			factorize_Symmetric_Non_Leaf(j, k);
		}
		qr_Level(j-1);
	}

	if(nLevels>0){
		small_Cholesky(0,0);
	}

	std::cout << "\nEnd Symmetric factorize...\n";
}

/******************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining QR decompositions of factors of low-rank approximations */
/* of the off-diagonal blocks at a certain level. It is used to ensure the U is unitary in (I+UKU'). 	          */
/******************************************************************************************************************/

void HODLR_Tree::qr_Level(int level){
    #pragma omp parallel for
	for(int k=0; k<nodesInLevel[level]; ++k){
		qr(level, k);
	}
}

/****************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining QR decomposition of a node at a level. 	*/
/****************************************************************************************************/

void HODLR_Tree::qr(int j, int k){

	int min0 = std::min(tree[j][k]->Qfactor[0].rows(), tree[j][k]->Qfactor[0].cols());
	int min1 = std::min(tree[j][k]->Qfactor[1].rows(), tree[j][k]->Qfactor[1].cols());

	Eigen::HouseholderQR<Eigen::MatrixXd> qr(tree[j][k]->Qfactor[0]);
	tree[j][k]->Qfactor[0] = qr.householderQ()*(Eigen::MatrixXd::Identity(tree[j][k]->Qfactor[0].rows(), min0));
	tree[j][k]->K = qr.matrixQR().block(0,0,min0,tree[j][k]->Qfactor[0].cols()).triangularView<Eigen::Upper>()*tree[j][k]->K;

	Eigen::HouseholderQR<Eigen::MatrixXd> qr1(tree[j][k]->Qfactor[1]);
	tree[j][k]->Qfactor[1] = qr1.householderQ()*(Eigen::MatrixXd::Identity(tree[j][k]->Qfactor[1].rows(), min1));
	tree[j][k]->K *= qr1.matrixQR().block(0,0,min1,tree[j][k]->Qfactor[1].cols()).triangularView<Eigen::Upper>().transpose();

}

/****************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining symmetric factorisations of leaf nodes	*/
/****************************************************************************************************/

void HODLR_Tree::factorize_Symmetric_Leaf(int k) {
	//	 std::cout << "\nStart Symmetric factorize Leaf: " << k << "\n";

	int parent	=	k;
	int child	=	k;
	int size	=	tree[nLevels][k]->nSize;
	int tstart, r;
    #pragma omp parallel for
	for (int l=nLevels-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[nLevels][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->sym_rank;
		tree[l][parent]->Qfactor[child].block(tstart,0,size,r) = solve_Symmetric_Factor_Leaf(k, tree[l][parent]->Qfactor[child].block(tstart,0,size,r));
	}
	//	 std::cout << "\nDone Symmetric factorize Leaf: " << k << "\n";
}

/********************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving Kx = b, where K is the symmetric factor of leaf node k. */
/********************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of K multiplied with b

Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Factor_Leaf(int k, Eigen::MatrixXd b) {

	return tree[nLevels][k]->llt.matrixL().solve(b);
}

/********************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is used for obtaining symmetric factorisations of non-leaf nodes	*/
/********************************************************************************************************/

void HODLR_Tree::factorize_Symmetric_Non_Leaf(int j, int k) {
	//	 std::cout << "\nStart Symmetric factorize; Level: " << j << " Node: " << k << "\n";
	small_Cholesky(j,k);
	int parent = k;
	int child	=	k;
	int size	=	tree[j][k]->nSize;
	int tstart, r;
    #pragma omp parallel for
	for (int l=j-1; l>=0; --l) {
		child	=	parent%2;
		parent	=	parent/2;
		tstart	=	tree[j][k]->nStart-tree[l][parent]->cStart[child];
		r		=	tree[l][parent]->sym_rank;
		tree[l][parent]->Qfactor[child].block(tstart,0,size,r)	=	solve_Symmetric_Factor_Non_Leaf(j, k, tree[l][parent]->Qfactor[child].block(tstart,0,size,r));
	}
	//	 std::cout << "\nDone factorize; Level: " << j << " Node: " << k << "\n";
}

/************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine is utilised in obtaining X in (I+UKU') = (I+UXU')(I+UXU')'.	*/
/************************************************************************************************/

void HODLR_Tree::small_Cholesky(int j, int k) {

	tree[j][k]->llt.compute(Eigen::MatrixXd::Identity(tree[j][k]->sym_rank, tree[j][k]->sym_rank) - tree[j][k]->K.transpose()*tree[j][k]->K);

}

/********************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving (I+UXU')x = b. The method uses Sherman-Morrison-Woodsbury formula to obtain x.	*/
/********************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	(I-U(inverse of (I+X))XU') multiplied by b
Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Factor_Non_Leaf(int j, int k, Eigen::MatrixXd b) {
	// std::cout << "\nStart Solve Symm non leaf; Level: " << j << "Node: " << k << "\n";
	int n0 = tree[j][k]->Q[0].rows();
	int n1 = tree[j][k]->Q[1].rows();
	//   std::cout<<"b.cols "<<b.cols()<<"\n";
	Eigen::MatrixXd tmp = tree[j][k]->Qfactor[1].transpose()*b.block(n0,0,n1,b.cols());
	b.block(n0,0,n1,b.cols()) -= tree[j][k]->Qfactor[1]*(tree[j][k]->llt.matrixL().solve((tree[j][k]->K.transpose()*(tree[j][k]->Qfactor[0].transpose()*b.block(0,0,n0,b.cols()))) - tmp) +tmp);
	return(b);
}

/************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains fast determinant of kernel matrix 	*/
/************************************************************************/

/************/
/*	OUTPUT	*/
/************/

///	determinant			-	Determinant of the kernel matrix

double HODLR_Tree::determinant(){
	double det = 0.0;
	for (int j=nLevels; j>=0; --j) {
		for (int k=0; k<nodesInLevel[j]; ++k) {
			for (int l=0; l<tree[j][k]->Kfactor.matrixLU().rows(); ++l) {
				det += log(fabs(tree[j][k]->Kfactor.matrixLU()(l,l)));
			}
		}
	}
	return(det);
}

/****************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains fast determinant of SPD kernel matrix whose symmetric factorization has been computed 	*/
/****************************************************************************************************************************/

/************/
/*	OUTPUT	*/
/************/

///	determinant			-	Determinant of the kernel matrix

double HODLR_Tree::symmetric_Determinant() {
	double det = 0.0;
	for (int j=nLevels; j>=0; --j) {
		for (int k=0; k<nodesInLevel[j]; ++k) {
			for (int l=0; l<tree[j][k]->llt.matrixL().rows(); ++l) {
				det += log(tree[j][k]->llt.matrixL()(l,l));
			}
		}
	}
	return 2*det;
}

/************************************************/
/*	PURPOSE OF EXISTENCE:	Solves for Wx = B	*/
/************************************************/

/************/
/*	INPUTS	*/
/************/

///	B			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of W multiplied by B

Eigen::MatrixXd HODLR_Tree::symmetric_Factor_Solve(Eigen::MatrixXd b1) {
	int start, size;
	Eigen::MatrixXd b = b1;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();
	Eigen::MatrixXd M1[nodesInLevel[nLevels]];
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		M1[k] = solve_Symmetric_Factor_Leaf(k, b.block(start, 0, size, r));
	}

	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r)	=  M1[k];
	}

	b=x;
	for (int j=nLevels-1; j>=0; --j) {
		Eigen::MatrixXd M2[nodesInLevel[j]];
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			M2[k] = solve_Symmetric_Factor_Non_Leaf(j, k, b.block(start, 0, size, r));
		}

		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	M2[k];
		}
		b=x;
	}

	return x;
}

/************************************************/
/*	PURPOSE OF EXISTENCE:	Solves for W'x = B	*/
/************************************************/

/************/
/*	INPUTS	*/
/************/

///	B			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of W' multiplied by B

Eigen::MatrixXd HODLR_Tree::symmetric_Factor_Transpose_Solve(Eigen::MatrixXd B) {
	int start, size;
	Eigen::MatrixXd b = B;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();
	Eigen::MatrixXd M1[nodesInLevel[nLevels]];
	for (int j=0; j<nLevels; ++j) {
		Eigen::MatrixXd M2[nodesInLevel[j]];
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			M2[k]	=	solve_Symmetric_Factor_Transpose_Non_Leaf(j, k, b.block(start, 0, size, r));
		}
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	M2[k];
		}
		b=x;
	}

    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		M1[k]	=	solve_Symmetric_Factor_Transpose_Leaf(k, b.block(start, 0, size, r));
	}

	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r) = M1[k];
	}

	return x;
}

/**********************************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Solves for the 	linear system Ax=b, where A is the kernel matrix whose symmetric factorization has been computed. */
/**********************************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of kernel matrix multiplied by input matrix

Eigen::MatrixXd HODLR_Tree::symmetric_Solve(Eigen::MatrixXd b) {
	return(symmetric_Factor_Transpose_Solve(symmetric_Factor_Solve(b)));
}

/************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving K'x = b, where K is the symmetric factor of leaf node k	*/
/************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of K' multiplied with b

Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Factor_Transpose_Leaf(int k, Eigen::MatrixXd b){
	return(tree[nLevels][k]->llt.matrixL().transpose().solve(b));
}

/************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for solving (I+UXU')'x = b	*/
/************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	x			-	Inverse of (I+UXU')' multiplied with b

Eigen::MatrixXd HODLR_Tree::solve_Symmetric_Factor_Transpose_Non_Leaf(int j, int k, Eigen::MatrixXd b){

	int n0 = tree[j][k]->Q[0].rows();
	int n1 = tree[j][k]->Q[1].rows();
	Eigen::MatrixXd xtmp = tree[j][k]->Qfactor[1].transpose()*b.block(n0,0,n1,b.cols());
	Eigen::MatrixXd ytmp = tree[j][k]->llt.matrixL().transpose().solve(xtmp);
	b.block(0, 0, n0, b.cols()) -= tree[j][k]->Qfactor[0]*(tree[j][k]->K*ytmp);
	b.block(n0, 0, n1, b.cols()) -= tree[j][k]->Qfactor[1]*(xtmp - ytmp);
	return(b);

}

/********************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains the product of W with an input matrix	*/
/********************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

///	B			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	W * B

Eigen::MatrixXd HODLR_Tree::symmetric_Factor_Product(Eigen::MatrixXd B){
	int start, size;
	Eigen::MatrixXd b = B;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();

	for (int j=0; j<nLevels; ++j) {
		Eigen::MatrixXd M2[nodesInLevel[j]];
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			M2[k]	=	mult_Symmetric_Factor_Non_Leaf(j, k, b.block(start, 0, size, r));
		}
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	M2[k];
		}
		b=x;
	}

	Eigen::MatrixXd M1[nodesInLevel[nLevels]];
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		M1[k] = tree[nLevels][k]->llt.matrixL()*b.block(start, 0, size, r);
	}

	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r)	= M1[k];
	}
	b = x;

	return(x);
}

/********************************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for obtaining the product of (I+UXU') with an input matrix	*/
/********************************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	(I+UXU') * b

Eigen::MatrixXd HODLR_Tree::mult_Symmetric_Factor_Non_Leaf(int j, int k, Eigen::MatrixXd b) {
	int n0 = tree[j][k]->Q[0].rows();
	int n1 = tree[j][k]->Q[1].rows();
	Eigen::MatrixXd tmp = tree[j][k]->Qfactor[1].transpose()*b.block(n0,0,n1,b.cols());
	b.block(n0,0,n1,b.cols()) += tree[j][k]->Qfactor[1]*((tree[j][k]->K.transpose()*tree[j][k]->Qfactor[0].transpose()*b.block(0,0,n0,b.cols())) + ((Eigen::MatrixXd)tree[j][k]->llt.matrixL() - Eigen::MatrixXd::Identity(tree[j][k]->sym_rank, tree[j][k]->sym_rank))*tmp);
	return(b);
}

/********************************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains the product of W' with an input matrix	*/
/********************************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

///	B			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	W' * B

Eigen::MatrixXd HODLR_Tree::symmetric_Factor_Transpose_Product(Eigen::MatrixXd B){
	int start, size;
	Eigen::MatrixXd b = B;
	Eigen::MatrixXd x	=	Eigen::MatrixXd::Zero(b.rows(),b.cols());
	int r	=	b.cols();


	Eigen::MatrixXd M1[nodesInLevel[nLevels]];
    #pragma omp parallel for
	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		M1[k] = tree[nLevels][k]->llt.matrixL().transpose()*b.block(start, 0, size, r);
	}

	for (int k=0; k<nodesInLevel[nLevels]; ++k) {
		start	=	tree[nLevels][k]->nStart;
		size	=	tree[nLevels][k]->nSize;
		x.block(start, 0, size, r)	= M1[k];
	}
	b = x;

	for (int j=nLevels-1; j>=0; --j) {
		Eigen::MatrixXd M2[nodesInLevel[j]];
        #pragma omp parallel for
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			M2[k]	=	mult_Symmetric_Factor_transpose_Non_Leaf(j, k, b.block(start, 0, size, r));
		}
		for (int k=0; k<nodesInLevel[j]; ++k) {
			start	=	tree[j][k]->nStart;
			size	=	tree[j][k]->nSize;
			x.block(start, 0, size, r)	=	M2[k];
		}
		b=x;
	}
	return(x);
}

/************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Routine for obtaining the product of (I+UXU')' with an input matrix	*/
/************************************************************************************************/

/************/
/*	INPUTS	*/
/************/

/// j           -   Level number
/// k           -   Node number
///	b			-	Input matrix

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	(I+UXU')' * b

Eigen::MatrixXd HODLR_Tree::mult_Symmetric_Factor_transpose_Non_Leaf(int j, int k, Eigen::MatrixXd b){
	int n0 = tree[j][k]->Q[0].rows();
	int n1 = tree[j][k]->Q[1].rows();
	Eigen::MatrixXd tmp = tree[j][k]->Qfactor[1].transpose()*b.block(n0,0,n1,b.cols());
	b.block(0,0,n0,b.cols()) += tree[j][k]->Qfactor[0]*(tree[j][k]->K*tmp);
	b.block(n0,0,n1,b.cols()) += tree[j][k]->Qfactor[1]*(((Eigen::MatrixXd)tree[j][k]->llt.matrixL().transpose() - Eigen::MatrixXd::Identity(tree[j][k]->sym_rank, tree[j][k]->sym_rank))*tmp);
	return(b);
}

/************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains the symmetric factor of the matrix	*/
/************************************************************************/

/************/
/*	OUTPUT	*/
/************/

///	matrix			-	symmetric factor of the matrix

Eigen::MatrixXd HODLR_Tree::build_Symmetric_Matrix_Factor() {

	if(nLevels == 0){
		return tree[0][0]->llt.matrixL();
	}

	Eigen::MatrixXd Wc;
	Eigen::MatrixXd Rc = Eigen::MatrixXd::Identity(N,N);
	for(int l=0;l<nLevels;++l){
		Wc = Eigen::MatrixXd::Zero(N,N);
        #pragma omp parallel for
		for(int k=0; k<nodesInLevel[l]; ++k){
			int r = tree[l][k]->sym_rank;
			int n0 = tree[l][k]->Q[0].rows();
			int n1 = tree[l][k]->Q[1].rows();
			Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n0+n1,n0+n1);
			T.block(n0,0,n1,n0) = tree[l][k]->Qfactor[1]*(tree[l][k]->K.transpose()*tree[l][k]->Qfactor[0].transpose());
			T.block(n0,n0,n1,n1) += tree[l][k]->Qfactor[1]*(((Eigen::MatrixXd)tree[l][k]->llt.matrixL() - Eigen::MatrixXd::Identity(r, r))*tree[l][k]->Qfactor[1].transpose());
			Wc.block(tree[l][k]->nStart, tree[l][k]->nStart, tree[l][k]->nSize, tree[l][k]->nSize) = T;
		}
		Rc = Wc*Rc;
	}

	Wc = Eigen::MatrixXd::Zero(N,N);

    #pragma omp parallel for
	for(int k=0; k<nodesInLevel[nLevels]; ++k){
		Wc.block(tree[nLevels-1][k/2]->cStart[k%2], tree[nLevels-1][k/2]->cStart[k%2], tree[nLevels-1][k/2]->cSize[k%2], tree[nLevels-1][k/2]->cSize[k%2]) = tree[nLevels][k]->llt.matrixL();
	}

	Rc = Wc*Rc;

	return Rc;
}

#endif /*__HODLR_Tree__*/
