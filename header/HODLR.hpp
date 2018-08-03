//version	:	1.1.1
//date		:	11-4-2018(dd-mm-yyyy)
//author	:	Vaishnavi G.

#ifndef _hodlr_HPP__
#define _hodlr_HPP__

#include <vector>
#include <Eigen/Dense>	
#include <Eigen/Sparse>
#include <Eigen/OrderingMethods>
#include <fstream>

#define EIGEN_DONT_PARALLELIZE
#include <iostream>
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include <time.h>   
using namespace Eigen;
using namespace std;

const double PI	=	3.1415926535897932384;
typedef Eigen::Triplet<double> T;

struct location {
	int row;
	int col;
};

class HODLR1DLine {
public:
	int lineNumber;
	
	HODLR1DLine () {
		lineNumber		=	-1;
	}
        
	double center;
	std::vector<double> chebNodes;
	std::vector<double> ChargeChebNodes;			// the interpolated charge locations at the chebynodes
};


class kernel {
public:
	double alpha;		//	Degree of homogeneity of the kernel.
	kernel() {};
	~kernel() {};
	virtual double getInteraction(const double r1, const double r2){
		return 0.0;
	};	//	Kernel entry generator
};


template <typename kerneltype>
class HODLR1DTree {
public:
	kerneltype* K;
	int nLevels;			//	Number of levels in the tree.
        int nChebNodes;		//=rank	//	Number of Chebyshev nodes along one direction.
	int N;				//	Number of particles.
	double L;			//	Semi-length of the simulation box.
	int rank;
	int sizeA;			//	size of the Extended Sparse Matrix
	Eigen::VectorXd b;	
	std::vector<int> nLinesPerLevel;//	Number of boxes at each level in the tree.
	std::vector<double> lineLength;	//	Box radius at each level in the tree assuming the box at the root is [-1,1]^2
	std::vector<double> shiftTo_nLevels;	//	As we go down the tree, the shift a block at a lower level undergoes when a new level comes in
	std::vector<double> newRowsPerLevel;	//	New rows that get added added to the matrix as we go down a level in the tree
	std::vector<std::vector<HODLR1DLine> > tree;	//	The tree storing all the information.
				
	//	Chebyshev nodes
	std::vector<double> standardChebNodes1D;
	std::vector<double> standardCharges1D;

	//	Operator
	Eigen::MatrixXd selfInteraction;		//	Needed only at the leaf level.

	Eigen::VectorXd extendedb;

	// Eigen Sparse Data Structure
	SparseMatrix<double, ColMajor> A;
	std::vector<T> tripletList;	
	
	HODLR1DTree(kerneltype* K, int nLevels, int nChebNodes, int N, const Eigen::VectorXd& b, Eigen::VectorXd& charges, double L) {
		this->K				=	K;
		this->b				=	b;
		this->L				=	L;
		this->nLevels			=	nLevels;
		this->nChebNodes		=	nChebNodes;//=rank
		L				=	L;
		this->rank			=	nChebNodes;
	        nLinesPerLevel.push_back(1);//level 0
		lineLength.push_back(L);//semi length
		shiftTo_nLevels.push_back(0); // data of level 0 not used;so a dummy value is placed
		newRowsPerLevel.push_back(2*rank); //dummy
		for (int k=1; k<=nLevels; ++k) {
			nLinesPerLevel.push_back(2*nLinesPerLevel[k-1]);
			lineLength.push_back(0.5*lineLength[k-1]);
			newRowsPerLevel.push_back(newRowsPerLevel[k-1]*2);
			if (k==1)
				shiftTo_nLevels.push_back(8*rank*(int(pow(2,nLevels-1))-1));
			else
				shiftTo_nLevels.push_back(shiftTo_nLevels[k-1]-newRowsPerLevel[k]);
			}
		this->N				=	N;
		this->sizeA			=	N+4*rank*(int(pow(2,nLevels))-1);
	}



        //	set_Standard_Cheb_Nodes
	void set_Standard_Cheb_Nodes() {
		for (int k=1; k<=nChebNodes; ++k) {
			standardChebNodes1D.push_back(-cos((k-0.5)/nChebNodes*PI));
		}
        }


	//	set_Standard_charges; (2*rank) charges in [-1,1]
	void set_Standard_Charges() {
		for (int k=1; k<=2*nChebNodes; ++k) {
			standardCharges1D.push_back(-cos((k-0.5)/(2*nChebNodes)*PI));
		}
        }


	//	shifted_scaled_cheb_nodes	
	void scalePoints(int NumPoints, std::vector<double> &oldPoints, double oldCenter, double oldRadius, std::vector<double> &newPoints, double newCenter, double newRadius) {
		for (int k=0; k<NumPoints; ++k) {
			double temp;
			temp	=	newRadius*(oldPoints[k]-oldCenter)/oldRadius+newCenter;
			newPoints.push_back(temp);
		}
	}


	//	get_ChebPoly
	double get_ChebPoly(double x, int n) {
		return cos(n*acos(x));
	}


	//	get_S
	double get_S(double x, double y, int n) {
		double S	=	0.5;
		for (int k=1; k<n; ++k) {
			S+=get_ChebPoly(x,k)*get_ChebPoly(y,k);
		}
		return 2.0/n*S;
	}


	// create tree data structure
	void createTree() {
		//	First create root and add to tree
		HODLR1DLine root;
		std::vector<HODLR1DLine> rootLevel;
		rootLevel.push_back(root);
		tree.push_back(rootLevel);

		for (int j=1; j<=nLevels; ++j) {
			std::vector<HODLR1DLine> level;
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				HODLR1DLine line;
				line.lineNumber		=	k;
				level.push_back(line);
			}
			tree.push_back(level);
		}
	}

	//assign centers of lines at all levels
	void assign_Center_Location() {
		int J, K;
		tree[0][0].center	=	0.0;		
		for (int j=0; j<nLevels; ++j) {
			J	=	j+1;
			double shift	=	0.5*lineLength[j];
			for (int k=0; k<nLinesPerLevel[j]; ++k) {
				K	=	2*k;
				tree[J][K].center	=	tree[j][k].center-shift;
				tree[J][K+1].center	=	tree[j][k].center+shift;
				
			}
		}
	}


	// assign charges at the leaf level; 2*rank charges per leaf
	void assign_Leaf_Charges() {
		for (int l=0; l<nLinesPerLevel[nLevels]; ++l) {//line
			scalePoints(2*nChebNodes, standardCharges1D, 0.0, 1.0, tree[nLevels][l].ChargeChebNodes, tree[nLevels][l].center,lineLength[nLevels]);
		}
	}
	

	//	assign ChebNodes at all non-zero levels and all lines
	void assign_ChebNodes() {
		for (int k=1; k<=nLevels; ++k) {//level
			for (int l=0; l<nLinesPerLevel[k]; ++l) {//line
				scalePoints(nChebNodes, standardChebNodes1D, 0.0, 1.0, tree[k][l].chebNodes, tree[k][l].center,lineLength[k]);
			}
		}	
	}


	// evalaute 'SelfInteractionOperator' at leaf level which is translation invariant
	void SelfInteractionOperator(Eigen::MatrixXd& M) {
		M	=	Eigen::MatrixXd(2*nChebNodes,2*nChebNodes);
		std::vector<double> scaledChargeChebNodes;
		scalePoints(2*nChebNodes, standardCharges1D, 0.0, 1.0, scaledChargeChebNodes, 0.0,lineLength[nLevels]);
		for (int i=0; i<2*nChebNodes; ++i) {
			for (int j=0; j<2*nChebNodes; ++j) {
				M(i,j)	=	K->getInteraction(scaledChargeChebNodes[i], scaledChargeChebNodes[j]);
			}
		}
	}

	

	// assemble the matrix in extended sparse format
	void assembleMatrix(int LowRankChoice) {
		int NNZ	=	5*N*N/int(pow(2,nLevels))+N*N*nLevels/int(pow(2,nLevels-1))+N*N; //no. of non-zero elements in the extended sparse matrix
		tripletList.reserve(NNZ); //reserving space for the NNZ elements 
		if (LowRankChoice == 1) {
			assembleUKV_ChebyshevInterpolation(); //assemble the elements corresponding to the U,V matrices at all levels, at their respective positions in the extended sparse matrix
		}
		else {
			assembleUKV_ACA();
		}
		assembleI(); //assemble -I
		assembleK(); //assemble self interaction at leaf level
	}


	// U: L2L; V: M2M
	void assembleUKV_ChebyshevInterpolation() {
		for (int k=1; k<=nLevels; ++k) {//level
			struct location offset; //offset location of U blocks (when there is a shift in blocks due to addition of higher level, this is not the global location)
			offset.row	=	0;
			offset.col	=	N;
			for (int l=0; l<nLinesPerLevel[k]; ++l) {//line
				std::vector<double> ChebyCharges; //assemble all charges that are present at a line in that level
				std::vector<double> scaledChebyCharges;
				for (int h=0; h<pow(2,nLevels-k); ++h) {
					for (int c=0; c< 2*nChebNodes; ++c) {
						ChebyCharges.push_back(tree[nLevels][l*pow(2,nLevels-k)+h].ChargeChebNodes[c]);
					}
				}
				scalePoints(2*nChebNodes*pow(2,nLevels-k), ChebyCharges, tree[k][l].center, lineLength[k],scaledChebyCharges, 0.0 ,1.0);
				for (int i=0; i<2*nChebNodes*pow(2,nLevels-k); ++i) {//number of rows in U=pow(2,nLevels-k)*nChebNodes
					for (int j=0; j<nChebNodes; ++j) {//the charges present on the line--j
						double U_ij	=	get_S(standardChebNodes1D[j], scaledChebyCharges[i], nChebNodes);
						tripletList.push_back(T(offset.row+i, shiftTo_nLevels[k]+offset.col+j, U_ij));
						tripletList.push_back(T(shiftTo_nLevels[k]+offset.col+j, offset.row+i, U_ij)); //using the fact that the matrix is symmetric
					}
				}
				offset.row	+=	2*nChebNodes*pow(2,nLevels-k);
				offset.col	+=	rank;
			}	
		}
		assembleK12_ChebyshevInterpolation(); //assemble low rank interaction among siblings
	}



	//K12: M2L
	void assembleK12_ChebyshevInterpolation() {
		for (int k=1; k<=nLevels; ++k) {//level
			struct location offset;
			offset.row	=	N+int(pow(2,k))*rank;//offset location of M2L blocks (when there is a shift in blocks due to addition of higher level, this is not the global location)
			offset.col	=	N+int(pow(2,k))*rank + rank;
			for (int l=0; l<nLinesPerLevel[k]; ++l) {//line
				int lother;
				if(l%2 == 0) {
					lother = l+1;
				}
				else {
					lother = l-1;
				}
				for (int i=0; i<nChebNodes; ++i) {//chebynodes--i
					for (int j=0; j<=i; ++j) {//the charges present on the line--j
						double K12_ij	=	K->getInteraction(tree[k][l].chebNodes[i], tree[k][lother].chebNodes[j]);
						tripletList.push_back(T(shiftTo_nLevels[k]+offset.row+i, shiftTo_nLevels[k]+offset.col+j, K12_ij));
						if(i != j) { //	using the fact that the matrix is symmetric
							tripletList.push_back(T(shiftTo_nLevels[k]+offset.col+j, shiftTo_nLevels[k]+offset.row+i, K12_ij));
						}
					}
				}
				offset.row	=	offset.row+rank;
				if (l%2 == 0) {
					offset.col	=	offset.col-rank;
				}
				else {
					offset.col	=	offset.col+3.0*rank;
				}
			}
		}
	}


		// ACA
	// ChebyCharges1: Line 1 charge locations (N1)
	// ChebyCharges2: Line 2 charge locations (N2)
	//when the matrix to be approximated is from a kernel which takes ChebyCharges1, ChebyCharges2 as the inputs
	void ACA(std::vector<double> ChebyCharges1, std::vector<double> ChebyCharges2, Eigen::MatrixXd& U, Eigen::MatrixXd& V, Eigen::VectorXd& K12) {
		// K12 is a diagonal matrix, so we are storing it as a vector
		int N1	=	ChebyCharges1.size();
		int N2	=	ChebyCharges2.size();
		U	=	Eigen::MatrixXd(N1,nChebNodes);
		V	=	Eigen::MatrixXd(nChebNodes,N2);
		K12	=	Eigen::VectorXd(nChebNodes);
		std::ptrdiff_t col_index;
		std::ptrdiff_t row_index;
		Eigen::RowVectorXd row;
		Eigen::VectorXd col;
		Eigen::RowVectorXd v;
		Eigen::VectorXd u;
		Eigen::RowVectorXd row_dumb;
		Eigen::VectorXd col_dumb;		
		// the matrix to be approximated has dimensions N1xN2
			
		// initialisation
		int k=0;
		Eigen::MatrixXd R	=	Eigen::MatrixXd::Zero(N1,N2);
		std::vector<double> row_list;
		std::vector<double> col_list;
		row_index=0;
		row_list.push_back(row_index);
		for (int g=0; g<N2; ++g) {
			R(row_index,g)	=	K->getInteraction(ChebyCharges1[row_index], ChebyCharges2[g]);
		}
		//R.row(i)=row;
		row	=	R.row(row_index).cwiseAbs();
		row.maxCoeff(&col_index);
		col_list.push_back(int(col_index));
		v=row/R(row_index,int(col_index));
		for (int g=0; g<N1; ++g) {
			R(g,col_index)	=	K->getInteraction(ChebyCharges1[g], ChebyCharges2[col_index]);
		}
		u	=	R.col(col_index);
		U.col(k)	=	u;
		V.row(k)	=	v;
		K12(k)		=	1.0;
		double normS	=	0.0;
		col_dumb	=	R.col(col_index);
		for (int e=0; e<row_list.size(); ++e) {
			col_dumb(row_list[e])	=	0.0;
		}		
		col	=	col_dumb.cwiseAbs();
		col.maxCoeff(&row_index);
		row_list.push_back(int(row_index));
		// iterations
		//while (u.norm()*v.norm() >= epsilon*normS) {
		while (k < nChebNodes-1) {
			++k;
			for (int g=0; g<N2; ++g) {
				R(row_index,g)	=	K->getInteraction(ChebyCharges1[row_index], ChebyCharges2[g]);
			}
			for (int l=0; l<=k-1; ++l){
				R.row(row_index)	-=	U(row_index,l)*V.row(l);
			}
			row_dumb	=	R.row(row_index);
			for (int e=0; e<col_list.size(); ++e) {
				row_dumb(col_list[e])	=	0.0;
			}		
			row	=	row_dumb.cwiseAbs();
			row.maxCoeff(&col_index);
			col_list.push_back(int(col_index));
			v=R.row(row_index)/R(row_index,col_index);
			for (int g=0; g<N1; ++g) {
				R(g,col_index)	=	K->getInteraction(ChebyCharges1[g], ChebyCharges2[col_index]);
			}
			for (int l=0; l<=k-1; ++l){
				R.col(col_index)	-=	V(l,col_index)*U.col(l);
			}
	
			u	=	R.col(col_index);
			U.col(k)	=	u;
			V.row(k)	=	v;
			K12(k)		=	1.0;
			/*normS		=	pow(normS,2)+pow(u.norm()*v.norm(),2);
			for (int l=0; l<=k-1; ++l){
				normS	+=	2*fabs(U.col(l).dot(u) * V.row(l).dot(v));
			}
			normS	=	sqrt(normS);*/
			col_dumb	=	R.col(col_index);
			for (int e=0; e<row_list.size(); ++e) {
				col_dumb(row_list[e])	=	0.0;
			}		
			col	=	col_dumb.cwiseAbs();
			col.maxCoeff(&row_index);
			row_list.push_back(int(row_index));
			
		} // end while
	} // end ACA


	// U: L2L; V: M2M
	// using ACA
	void assembleUKV_ACA() {
		for (int k=1; k<=nLevels; ++k) {//level
			struct location offset; //offset location of U blocks (when there is a shift in blocks due to addition of higher level, this is not the global location)
			offset.row	=	0;
			offset.col	=	N;
			for (int l=0; l<nLinesPerLevel[k]; l=l+2) {//line
				std::vector<double> ChebyCharges1; //assemble all charges that are present at a line in that level
				std::vector<double> ChebyCharges2;
				for (int h=0; h<pow(2,nLevels-k); ++h) {
					for (int c=0; c< 2*nChebNodes; ++c) {
						ChebyCharges1.push_back(tree[nLevels][l*pow(2,nLevels-k)+h].ChargeChebNodes[c]);
						ChebyCharges2.push_back(tree[nLevels][(l+1)*pow(2,nLevels-k)+h].ChargeChebNodes[c]);
					}
				}
				Eigen::MatrixXd U;
				Eigen::MatrixXd V;
				Eigen::VectorXd K12;
				ACA(ChebyCharges1, ChebyCharges2, U, V, K12);//when the matrix to be approximated is from a kernel which takes ChebyCharges1, ChebyCharges2 as the inputs
				// Putting U into the extended sparse matrix
				for (int i=0; i<2*nChebNodes*pow(2,nLevels-k); ++i) {//number of rows in U=pow(2,nLevels-k)*nChebNodes
					for (int j=0; j<nChebNodes; ++j) {//the charges present on the line--j
						tripletList.push_back(T(offset.row+i, shiftTo_nLevels[k]+offset.col+j, U(i,j))); //U
						tripletList.push_back(T(shiftTo_nLevels[k]+offset.col+j, offset.row+i, U(i,j))); //V=U' //using the fact that the matrix is symmetric //U=V
					}
				}
				offset.row	+=	2*nChebNodes*pow(2,nLevels-k);
				offset.col	+=	rank;
				// Putting V into the extended sparse matrix
				for (int i=0; i<2*nChebNodes*pow(2,nLevels-k); ++i) {//number of rows in U=pow(2,nLevels-k)*nChebNodes
					for (int j=0; j<nChebNodes; ++j) {//the charges present on the line--j
						tripletList.push_back(T(shiftTo_nLevels[k]+offset.col+j, offset.row+i, V(j,i)));//V
						tripletList.push_back(T(offset.row+i, shiftTo_nLevels[k]+offset.col+j, V(j,i)));//U=V' //using the fact that the matrix is symmetric //U=V
					}
				}
				
				offset.row	+=	2*nChebNodes*pow(2,nLevels-k);
				offset.col	+=	rank;
				// Assemble K12
				assembleK12_ACA(k, l, K12);
			}	
		}
	}


	
		// M2L; Using ACA
		void assembleK12_ACA(int k, int l, Eigen::VectorXd& K12) {
			struct location offset;
			offset.row	=	N+int(pow(2,k))*rank;//offset location of M2L blocks (when there is a shift in blocks due to addition of higher level, this is not the global location)
			offset.col	=	N+int(pow(2,k))*rank;

			//l=0; K12
			offset.row	=	offset.row+(l/2)*(2*rank);
			offset.col	=	offset.col+(l/2)*(2*rank)+rank;
			for (int i=0; i<nChebNodes; ++i) {
				tripletList.push_back(T(shiftTo_nLevels[k]+offset.row+i, shiftTo_nLevels[k]+offset.col+i, K12(i)));
			}
			//K12=K21 (R symmetry)
			//l=1; K21
			offset.row	=	offset.row+rank;
			offset.col	=	offset.col-rank;
			for (int i=0; i<nChebNodes; ++i) {
				tripletList.push_back(T(shiftTo_nLevels[k]+offset.row+i, shiftTo_nLevels[k]+offset.col+i, K12(i)));
			}
		}



	void assembleI() {
		for (int k=1; k<=nLevels; ++k) {//level
			struct location offset;
			offset.row	=	N; //offset location of -I blocks (when there is a shift in blocks due to addition of higher level, this is not the global location)
			offset.col	=	N+pow(2,k)*rank;
			for (int l=0; l<nLinesPerLevel[k]; ++l) {//line
				for (int i=0; i<nChebNodes; ++i) {//chebynodes--i
					tripletList.push_back(T(shiftTo_nLevels[k]+offset.row+i, shiftTo_nLevels[k]+offset.col+i, -1.0));
					tripletList.push_back(T(shiftTo_nLevels[k]+offset.col+i, shiftTo_nLevels[k]+offset.row+i, -1.0));
				}
				offset.row	+=	rank; 
				offset.col	+=	rank;
			}		
		}
	}
	



	void assembleK() {
		struct location offset;
		offset.row	=	0;//offset location of U blocks // there is no shift occuring for the K blocks; so this is the global block location
		offset.col	=	0;
		Eigen::MatrixXd SelfInteraction;
		SelfInteractionOperator(SelfInteraction);
		for (int l=0; l<nLinesPerLevel[nLevels]; ++l) {//line
			for (int i=0; i<2*nChebNodes; ++i) {
				for (int j=0; j<=i; ++j) {
					double K_ij	=	SelfInteraction(i,j);
					tripletList.push_back(T(offset.row+i, offset.col+j, K_ij));
					if (i != j) { //	using the fact that the matrix is symmetric
						tripletList.push_back(T(offset.col+j, offset.row+i, K_ij));
					}
				}
			}
			offset.row	+=	2*nChebNodes;
			offset.col	+=	2*nChebNodes;
		}
	}



	void solve() { //	Using Eigen Sparse Solver
		A.resize(sizeA,sizeA);
		A.setFromTriplets(tripletList.begin(), tripletList.end());
		//cout << "A: " << endl << A << endl;
		SparseLU<SparseMatrix<double, ColMajor>, COLAMDOrdering<int> >   solver;
		//SparseLU<SparseMatrix<double, ColMajor>, NaturalOrdering<int> >   solver;
		solver.compute(A);
		extendedb	=	Eigen::VectorXd::Zero(sizeA);
		extendedb.head(N) = b;
		K->charges = solver.solve(extendedb);
	}

	
	double perform_Error_Check(int lineNumber) { //Performing Error Check for a particular Line
		Eigen::MatrixXd Aline(2*nChebNodes,N);
		for (int j=0; j<2*nChebNodes; ++j) {
			int i = 0;
			for (int l=0; l<pow(2,nLevels); ++l) {
				for (int c=0; c<2*nChebNodes; ++c) {
					Aline(j,i) = K->getInteraction(tree[nLevels][l].ChargeChebNodes[c],tree[nLevels][lineNumber].ChargeChebNodes[j]);
					++i;
				}
			}
		}
		
		Eigen::VectorXd bline(2*nChebNodes);
		for (int h=0; h<2*nChebNodes; ++h) {
			bline(h) = b(lineNumber*2*nChebNodes+h);
		}
		Eigen::VectorXd calculatedb(2*nChebNodes);
		Eigen::VectorXd abs_bline(2*nChebNodes);
		calculatedb	=	Aline*K->charges.head(N);
		Eigen::VectorXd error(2*nChebNodes);
		for (int k=0; k<2*nChebNodes; ++k) {
			error(k)	=	fabs((bline-calculatedb)(k));
			abs_bline(k)		=	fabs(bline(k));
		}
		return error.maxCoeff()/abs_bline.maxCoeff();
	}
};
#endif
