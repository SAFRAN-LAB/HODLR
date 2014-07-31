/*!
 \class HODLR_Node
 
 \brief This class is for the node of a HODLR tree, i.e., a diagonal sub-matrix at the appropriate level in the tree.

 \note
 
 \author $Sivaram Ambikasaran$
 
 \version
 
 \date $November 8th, 2013$
 
 Contact: siva.1985@gmail.com
 */

#ifndef __HODLR_NODE_HPP__
#define __HODLR_NODE_HPP__

#include <vector>
#include <Eigen/Dense>

using std::vector;
using namespace Eigen;

template <typename MatrixType>
class HODLR_Node {
private:
	MatrixType* kernel;

public:
	HODLR_Node* parent;
	HODLR_Node* child[2];

        /**	Variables of interest for both leaf and non-leaf;*/
	int levelNumber;           ///<	Level number of the node;
	int nodeNumber;            ///<	Node number is either 0 or 1;
        int levelBasedNodeNumber;       ///<    Level based node number;
	int nStart;                ///<	nStart is the starting index of node;
	int nSize;                 ///<	nSize is the size of the node;
	MatrixXd K;                     ///<	At leaf level, stores the self interaction; At non-leaf stores the matrix [I, V1inverse*U1inverse; V0inverse*U0inverse, I];
	FullPivLU<MatrixXd> Kinverse;	///<	Stores Factorization of K;
	double determinant;             ///<	Stores K.determinant();

        /**	Variables of interest for non-leaf;*/
	int nRank[2];	///<	nRank[0] rank of K01; nRank[1] rank of K10;
	MatrixXd U[2];		///<	Column basis of low-rank interaction of children at non-leaf;
	MatrixXd V[2];		///<	Row basis of low-rank interaction ofchildren at non-leaf;
	MatrixXd Uinverse[2];	///<	Column basis of low-rank interaction of children of inverse at non-leaf.
	MatrixXd Vinverse[2];	///<	Row basis of low-rank interaction of children of inverse at non-leaf.

        /**	Variables of interest at leaf;*/
	bool isLeaf;		///<	If node is a leaf, it takes the value TRUE;

        /*!
         Constructor for the class.
         */
	HODLR_Node(MatrixType* kernel, int levelNumber, int nodeNumber, int levelBasedNodeNumber, int nStart, int nSize) {
		this->kernel                    =	kernel;
		this->levelNumber               =	levelNumber;
                this->levelBasedNodeNumber      =       levelBasedNodeNumber;
		this->nodeNumber                =	nodeNumber;
		this->nStart                    =	nStart;
		this->nSize                     =	nSize;
		this->parent                    =	NULL;
		this->child[0]                  =	NULL;
		this->child[1]                  =	NULL;
	};

        /*!
         Assemble the relevant matrices.
         */
	void assemble_Matrices(double lowRankTolerance, VectorXd& diagonal, char s) {
		if (isLeaf	==	true) {
			kernel->get_Matrix(nStart, nStart, nSize, nSize, K);
			for (int k=0; k<nSize; ++k) {
				K(k,k)	=	diagonal(nStart+k);
			}
		}
		else if (isLeaf	==	false) {
			partial_Piv_LU(child[0]->nStart, child[1]->nStart, child[0]->nSize, child[1]->nSize, lowRankTolerance, nRank[0], U[0], V[1]);
            if (s == 's') {
                V[0]    =   U[0].transpose();
                U[1]    =   V[1].transpose();
                nRank[1]=   nRank[0];
            }
            else {
			partial_Piv_LU(child[1]->nStart, child[0]->nStart, child[1]->nSize, child[0]->nSize, lowRankTolerance, nRank[1], U[1], V[0]);
            }
//                        ranks[levelNumber][levelBasedNodeNumber][0]     =       nRank[0];
//                        ranks[levelNumber][levelBasedNodeNumber][1]     =       nRank[1];
		}
	};

        /*!
         Matrix matrix product.
         */
	void matrix_Matrix_Product(MatrixXd& x, MatrixXd& b) {
		int n	=	x.cols();

		if (isLeaf	==	true) {
			b.block(nStart, 0, nSize, n)	=	b.block(nStart, 0, nSize, n)	+	K*x.block(nStart, 0, nSize, n);
		}
		else if (isLeaf	==	false) {
			int n0	=	child[0]->nStart;
			int n1	=	child[1]->nStart;
			int m0	=	child[0]->nSize;
			int m1	=	child[1]->nSize;
			b.block(n0, 0, m0, n)	=	b.block(n0, 0, m0, n)	+	U[0]*(V[1]*x.block(n1, 0, m1, n));
			b.block(n1, 0, m1, n)	=	b.block(n1, 0, m1, n)	+	U[1]*(V[0]*x.block(n0, 0, m0, n));
		}
	};

        /*!
         Set UV Inversion.
         */
	void set_UV_Inversion() {
		for (int k=0; k<2; ++k) {
			Uinverse[k]	=	U[k];
			Vinverse[k]	=	V[k];
		}
	};

	void compute_K(char s) {
		if (isLeaf	==	false) {
			int m0	=	V[0].rows();
			int m1	=	V[1].rows();
			K	=	MatrixXd::Identity(m0+m1, m0+m1);

			K.block(0, m1, m1, m0)	=	Vinverse[1]*Uinverse[1];
            if (s=='s') {
//                K.block(m1, 0, m0, m1)  =   K.block(0, m1, m1, m0).transpose();
                K.block(m1, 0, m0, m1)	=	Vinverse[0]*Uinverse[0];
            }
            else {
                K.block(m1, 0, m0, m1)	=	Vinverse[0]*Uinverse[0];
            }
		}
	};

        /*!
         Computes the inverse.
         */
	void compute_Inverse() {
		Kinverse.compute(K);
	};

        /*!
         Applies the inverse.
         */
	void apply_Inverse(MatrixXd& matrix, int mStart) {
		int n	=	matrix.cols();
		int start	=	nStart-mStart;
		if (isLeaf	==	true) {
			matrix.block(start, 0, nSize, n)	=	Kinverse.solve(matrix.block(start, 0, nSize, n));
		}
		else if (isLeaf	==	false) {
			//	Computes temp		=	Vinverse*matrix

			MatrixXd temp(nRank[0]+nRank[1], n);

			temp.block(0, 0, nRank[0] , n)		=	Vinverse[1]*matrix.block(start+child[0]->nSize, 0 , child[1]->nSize, n);

			temp.block(nRank[0], 0, nRank[1] , n)	=	Vinverse[0]*matrix.block(start, 0 , child[0]->nSize, n);

			//	Computes tempSolve	=	Kinverse\temp

			MatrixXd tempSolve	=	Kinverse.solve(temp);

			//	Computes matrix		=	matrix-Uinverse*tempSolve

			matrix.block(start, 0, child[0]->nSize, n)			=	matrix.block(start, 0, child[0]->nSize, n)	-	Uinverse[0]*tempSolve.block(0, 0, nRank[0], n);
			matrix.block(start + child[0]->nSize, 0, child[1]->nSize, n)	=	matrix.block(start + child[0]->nSize, 0, child[1]->nSize, n)	-	Uinverse[1]*tempSolve.block(nRank[0], 0, nRank[1], n);
		}
	};

        /*!
         Computes the determinant of the matrix.
         */
	void compute_Determinant() {
                if (Kinverse.rows()>0) {        //      Check needed when the matrix is predominantly diagonal.
                        MatrixXd LU     =       Kinverse.matrixLU();
                        determinant     =       log(fabs(LU(0,0)));
                        for (int k=1; k<Kinverse.rows(); ++k) {
                                determinant+=log(fabs(LU(k,k)));
                        }
                        //              Previous version which had some underflow.
                        //              determinant	=	log(fabs(K.determinant()));
                }
	};

        /*!
         Destructor.
         */
	~HODLR_Node() {
		delete child[0];
		delete child[1];
		child[0]	=	NULL;
		child[1]	=	NULL;
	};

        /*!
         \brief Partial pivoted LU to construct low-rank.


         */
	void partial_Piv_LU(const int start_Row, const int start_Col, const int n_Rows, const int n_Cols, const double tolerance, int& computed_Rank, MatrixXd& U, MatrixXd& V) {

	/********************************/
	/*	PURPOSE OF EXISTENCE	*/
	/********************************/

	/*!
         Obtains the low-rank decomposition of the matrix to a desired tolerance using the partial pivoting LU algorithm, i.e., given a sub-matrix 'A' and tolerance 'epsilon', computes matrices 'U' and 'V' such that ||A-UV||_F < epsilon. The norm is Frobenius norm.
         */

	/************************/
	/*	INPUTS          */
	/************************/

	///	start_Row	-	Starting row of the sub-matrix.
	///	start_Col	-	Starting column of the sub-matrix.
	///	n_Rows		-	Number of rows of the sub-matrix.
	///	n_Cols		-	Number of columns of the sub-matrix.
	///	tolerance	-	Tolerance of low-rank approximation.

	/************************/
	/*	OUTPUTS		*/
	/************************/

	///	computed_Rank	-	Rank obtained for the given tolerance.
	///	U		-	Matrix forming the column basis.
	///	V		-	Matrix forming the row basis.

		/// If the matrix is small enough, do not do anything
		int tolerable_Rank =   5;
		if (n_Cols <= tolerable_Rank){
			kernel->get_Matrix(start_Row, start_Col, n_Rows, n_Cols, U);
			V               =   MatrixXd::Identity(n_Cols, n_Cols);
			computed_Rank   =   n_Cols;
			return;
		}
		else if (n_Rows <= tolerable_Rank){
			U               =   MatrixXd::Identity(n_Rows, n_Rows);
			kernel->get_Matrix(start_Row, start_Col, n_Rows, n_Cols, V);
			computed_Rank   =   n_Rows;
			return;
		}

		vector<int> rowIndex;	///	This stores the row indices, which have already been used.
		vector<int> colIndex;	///	This stores the column indices, which have already been used.
		vector<VectorXd> u;	///	Stores the column basis.
		vector<VectorXd> v;	///	Stores the row basis.

		srand (time(NULL));
		double max, Gamma, unused_max;

		/*  INITIALIZATION  */

		/// Initialize the matrix norm and the the first row index
		double matrix_Norm  =   0;
		rowIndex.push_back(0);

		int pivot;

		computed_Rank   =   0;

		VectorXd a, row, col;

		double row_Squared_Norm, row_Norm, col_Squared_Norm, col_Norm;

		/// Repeat till the desired tolerance is obtained
		do {
			/// Generation of the row
			kernel->get_Matrix_Row(start_Col, n_Cols, start_Row+rowIndex.back(), a);
			/// Row of the residuum and the pivot column
			row =   a;
			for (int l=0; l<computed_Rank; ++l) {
				row =   row-u[l](rowIndex.back())*v[l];
			}

			pivot   =   kernel->max_Abs_Vector(row, colIndex, max);

			int max_tries  =   100;
			int count      =   0;
			int count1     =   0;

			/// This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows.
			while (fabs(max)<tolerance && count < max_tries) {
				int new_rowIndex;
				rowIndex.pop_back();
				do {
					new_rowIndex   =   rand()%n_Rows;
					++count1;
				} while (find(rowIndex.begin(),rowIndex.end(),new_rowIndex)!=rowIndex.end() && count1 < max_tries);
				count1  =   0;
				rowIndex.push_back(new_rowIndex);

				/// Generation of the row
				kernel->get_Matrix_Row(start_Col, n_Cols, start_Row+rowIndex.back(), a);

				/// Row of the residuum and the pivot column
				row =   a;
				for (int l=0; l<computed_Rank; ++l) {
					row =   row-u[l](rowIndex.back())*v[l];
				}
				pivot   =   kernel->max_Abs_Vector(row, colIndex, max);
				++count;
			}

			if (count == max_tries) break;

			count = 0;

			colIndex.push_back(pivot);

			/// Normalizing constant
			Gamma   =   1.0/max;

			/// Generation of the column
			kernel->get_Matrix_Col(start_Row, n_Rows, start_Col+colIndex.back(), a);

			/// Column of the residuum and the pivot row
			col =   a;
			for (int l=0; l<computed_Rank; ++l) {
				col =   col-v[l](colIndex.back())*u[l];
			}
			pivot   =   kernel->max_Abs_Vector(col, rowIndex, unused_max);

			/// This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
			while (fabs(max)<tolerance && count < max_tries) {
				colIndex.pop_back();
				int new_colIndex;
				do {
					new_colIndex   =   rand()%n_Cols;
				} while (find(colIndex.begin(),colIndex.end(),new_colIndex)!=colIndex.end() && count1 < max_tries);
				count1  =   0;
				colIndex.push_back(new_colIndex);

				/// Generation of the column
				kernel->get_Matrix_Col(start_Row, n_Rows, start_Col+colIndex.back(), a);

				/// Column of the residuum and the pivot row
				col =   a;
				for (int l=0; l<computed_Rank; ++l) {
					col =   col-u[l](colIndex.back())*v[l];
				}
				pivot   =   kernel->max_Abs_Vector(col, rowIndex, unused_max);
				++count;
			}

			if (count == max_tries) break;

			count = 0;

			rowIndex.push_back(pivot);

			/// New vectors
			u.push_back(Gamma*col);
			v.push_back(row);

			/// New approximation of matrix norm
			row_Squared_Norm    =   row.squaredNorm();
			row_Norm            =   sqrt(row_Squared_Norm);

			col_Squared_Norm    =   col.squaredNorm();
			col_Norm            =   sqrt(col_Squared_Norm);

			matrix_Norm         =   matrix_Norm +   Gamma*Gamma*row_Squared_Norm*col_Squared_Norm;

			for (int j=0; j<computed_Rank; ++j) {
				matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
			}
			++computed_Rank;
		} while (row_Norm*col_Norm > fabs(max)*tolerance*matrix_Norm && computed_Rank <= fmin(n_Rows, n_Cols));

		/// If the computed_Rank is close to full-rank then return the trivial full-rank decomposition
		if (computed_Rank>=fmin(n_Rows, n_Cols)) {
			if (n_Rows < n_Cols) {
				U   =   MatrixXd::Identity(n_Rows,n_Rows);
				kernel->get_Matrix(start_Row, start_Col, n_Rows, n_Cols, V);
				computed_Rank   =   n_Rows;
				return;
			}
			else {
				kernel->get_Matrix(start_Row, start_Col, n_Rows, n_Cols, U);
				V   =   MatrixXd::Identity(n_Cols,n_Cols);
				computed_Rank   =   n_Cols;
				return;
			}
		}

		U   =   MatrixXd(n_Rows,computed_Rank);
		V   =   MatrixXd(computed_Rank,n_Cols);
		for (int j=0; j<computed_Rank; ++j) {
			U.col(j)    =   u[j];
			V.row(j)    =   v[j];
		}
	};

};

#endif /* defined(__HODLR_NODE_HPP__) */
