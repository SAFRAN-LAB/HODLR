#ifndef __HODLR_Matrix__
#define __HODLR_Matrix__

#include <Eigen/Dense>
#include <set>
#include <vector>

class HODLR_Matrix {
	friend class HODLR_Tree;
private:
	int N;
public:

	HODLR_Matrix(int N);

	virtual double get_Matrix_Entry(int j, int k) {
		return 0.0;
	}

	Eigen::VectorXd get_Row(int j, int nColStart, int nColEnd);
	Eigen::VectorXd get_Col(int k, int nRowStart, int nRowEnd);
	Eigen::MatrixXd get_Matrix(int j, int k, int nRows, int nCols);
	int max_Abs_Vector(const Eigen::VectorXd& v, const std::set<int>& allowed_Indices, double& max);
	void rook_Piv(int nRowStart, int nColStart, int nRows, int nCols, double tolerance, Eigen::MatrixXd& L, Eigen::MatrixXd& R, int& computedRank);
	void queen_Piv(int nRowStart, int nColStart, int nRows, int nCols, double tolerance, Eigen::MatrixXd& L, Eigen::MatrixXd& R, int& computedRank);
	~HODLR_Matrix();
};

/********************************************************/
/*	PURPOSE OF EXISTENCE:	Constructor for the class	*/
/********************************************************/

/************/
/*	INPUTS	*/
/************/

///	N	-	Number of rows/columns of the matrix

HODLR_Matrix::HODLR_Matrix(int N) {
	this->N	=	N;
}

/********************************************************/
/*	PURPOSE OF EXISTENCE:	Destructor for the class	*/
/********************************************************/

HODLR_Matrix::~HODLR_Matrix() {};

/********************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains a row of the matrix	*/
/********************************************************/

/************/
/*	INPUTS	*/
/************/

///	j			-	row index
///	nColStart	-	Starting column
///	nCols		-	Number of columns

/************/
/*	OUTPUT	*/
/************/

///	row			-	row of the matrix

Eigen::VectorXd HODLR_Matrix::get_Row(const int j, const int nColStart, const int nCols) {
	Eigen::VectorXd row(nCols);
	for (int k=0; k<nCols; ++k) {
		row(k)	=	get_Matrix_Entry(j,k+nColStart);
	}
	return row;
}

/************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains a column of the matrix	*/
/************************************************************/

/************/
/*	INPUTS	*/
/************/

///	k			-	column index
///	nRowStart	-	Starting row
///	nRows		-	Number of rows

/************/
/*	OUTPUT	*/
/************/

///	col			-	column of the matrix

Eigen::VectorXd HODLR_Matrix::get_Col(const int k, const int nRowStart, const int nRows) {
	Eigen::VectorXd col(nRows);
	for (int j=0; j<nRows; ++j) {
		col(j)	=	get_Matrix_Entry(j+nRowStart,k);
	}
	return col;
}

/****************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains a sub-matrix of the matrix	*/
/****************************************************************/

/************/
/*	INPUTS	*/
/************/

///	nRowStart	-	starting row
///	nColStart	-	Starting column
///	nRows		-	Number of rows
/// nCols		-	Number of columns

/************/
/*	OUTPUT	*/
/************/

///	mat			-	sub-matrix of the matrix

Eigen::MatrixXd HODLR_Matrix::get_Matrix(const int nRowStart, const int nColStart, const int nRows, const int nCols) {
	Eigen::MatrixXd mat(nRows, nCols);
	for (int j=0; j<nRows; ++j) {
		for (int k=0; k<nCols; ++k) {
			mat(j,k)	=	get_Matrix_Entry(j+nRowStart, k+nColStart);
		}
	}
	return mat;
}

/********************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains the index and value of the maximum entry of a vector 	*/
/********************************************************************************************/

/************/
/*	INPUTS	*/
/************/

///	v				-	vector
///	allowed_Indices	-	indices that need to be searched

/************/
/*	OUTPUT	*/
/************/

///	max		-	Value of the maximum
/// index	-	Index of the maximum

int HODLR_Matrix::max_Abs_Vector(const Eigen::VectorXd& v, const std::set<int>& allowed_Indices, double& max) {
	std::set<int>::iterator it	=	allowed_Indices.begin();
	int index	=	*allowed_Indices.begin();
	max				=	v(index);
	for (it=allowed_Indices.begin(); it!=allowed_Indices.end(); ++it) {
		if(fabs(v(*it))>fabs(max)) {
			index	=	*it;
			max		=	v(index);
		}
	}
	return index;
}

/****************************************************************************************************/
/*	PURPOSE OF EXISTENCE:	Obtains the low-rank decomposition of the matrix to a desired tolerance	*/
/* 	using rook pivoting, i.e., given a sub-matrix 'A' and tolerance 'epsilon', computes matrices	*/
/*	'L' and 'R' such that ||A-LR'||_F < epsilon. The norm is Frobenius norm.						*/
/****************************************************************************************************/

void HODLR_Matrix::rook_Piv(int nRowStart, int nColStart, int nRows, int nCols, double tolerance, Eigen::MatrixXd& L, Eigen::MatrixXd& R, int& computedRank) {
	std::vector<int> rowIndex;		///	This stores the row indices, which have already been used.
	std::vector<int> colIndex;		///	This stores the column indices, which have already been used.
	std::set<int> remainingRowIndex;/// Remaining row indicies
	std::set<int> remainingColIndex;/// Remaining row indicies
	std::vector<Eigen::VectorXd> u;	///	Stores the column basis.
	std::vector<Eigen::VectorXd> v;	///	Stores the row basis.

	for (int k=0; k<nRows; ++k) {
		remainingRowIndex.insert(k);
	}
	for (int k=0; k<nCols; ++k) {
		remainingColIndex.insert(k);
	}

	srand (time(NULL));
	double max, Gamma, unused_max;

	/*  INITIALIZATION  */

	/// Initialize the matrix norm and the the first row index
	double matrix_Norm  =   0;
	rowIndex.push_back(0);
	remainingRowIndex.erase(0);

	int pivot;

	computedRank   =   0;

	Eigen::VectorXd a, row, col;

	double row_Squared_Norm, row_Norm, col_Squared_Norm, col_Norm;

	int max_tries  =   10;

	int count;

	/// Repeat till the desired tolerance is obtained
	do {
		/// Generation of the row
		/// Row of the residuum and the pivot column
		// row =   A.row(rowIndex.back());
		row	=	get_Row(nRowStart+rowIndex.back(), nColStart, nCols);
		for (int l=0; l<computedRank; ++l) {
			row =   row-u[l](rowIndex.back())*v[l];
		}

		pivot   =   max_Abs_Vector(row, remainingColIndex, max);


		count	=   0;

		/// This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows upto some tolerance.
		while (fabs(max)<tolerance && count < max_tries && remainingColIndex.size() >0 && remainingRowIndex.size() >0) {
			rowIndex.pop_back();
			int new_rowIndex	=	*remainingRowIndex.begin();
			rowIndex.push_back(new_rowIndex);
			remainingRowIndex.erase(new_rowIndex);

			/// Generation of the row
			// a	=	A.row(new_rowIndex);
			a	=	get_Row(nRowStart+new_rowIndex, nColStart, nCols);
			/// Row of the residuum and the pivot column
			row =   a;
			for (int l=0; l<computedRank; ++l) {
				row =   row-u[l](rowIndex.back())*v[l];
			}
			pivot   =   max_Abs_Vector(row, remainingColIndex, max);
			++count;
		}

		if (count == max_tries || remainingColIndex.size() == 0 || remainingRowIndex.size() == 0) break;

		count = 0;

		colIndex.push_back(pivot);
		remainingColIndex.erase(pivot);

		/// Normalizing constant
		Gamma   =   1.0/max;

		/// Generation of the column
		// a	=	A.col(colIndex.back());
		a	=	get_Col(nColStart+colIndex.back(), nRowStart, nRows);
		/// Column of the residuum and the pivot row
		col =   a;
		for (int l=0; l<computedRank; ++l) {
			col =   col-v[l](colIndex.back())*u[l];
		}
		pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);

		/// This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
		while (fabs(max)<tolerance && count < max_tries && remainingColIndex.size() >0 && remainingRowIndex.size() >0) {
			colIndex.pop_back();
			int new_colIndex	=	*remainingColIndex.begin();
			colIndex.push_back(new_colIndex);
			remainingColIndex.erase(new_colIndex);

			/// Generation of the column
			// a	=	A.col(new_colIndex);
			a	=	get_Col(nColStart+new_colIndex, nRowStart, nRows);

			/// Column of the residuum and the pivot row
			col =   a;
			for (int l=0; l<computedRank; ++l) {
				col =   col-u[l](colIndex.back())*v[l];
			}
			pivot   =   max_Abs_Vector(col, remainingRowIndex, unused_max);
			++count;
			std::cout << count << "\n";
		}

		if (count == max_tries || remainingColIndex.size() == 0 || remainingRowIndex.size() == 0) break;

		count = 0;

		rowIndex.push_back(pivot);
		remainingRowIndex.erase(pivot);

		/// New vectors
		u.push_back(Gamma*col);
		v.push_back(row);

		/// New approximation of matrix norm
		row_Squared_Norm    =   row.squaredNorm();
		row_Norm            =   sqrt(row_Squared_Norm);

		col_Squared_Norm    =   col.squaredNorm();
		col_Norm            =   sqrt(col_Squared_Norm);

		matrix_Norm         =   matrix_Norm +   Gamma*Gamma*row_Squared_Norm*col_Squared_Norm;

		for (int j=0; j<computedRank; ++j) {
			matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
		}
		++computedRank;
	} while (computedRank*(nRows+nCols)*row_Norm*col_Norm > fabs(max)*tolerance*matrix_Norm && computedRank < fmin(nRows, nCols));

	/// If the computedRank is close to full-rank then return the trivial full-rank decomposition

	if (computedRank>=fmin(nRows, nCols)-1) {
		if (nRows < nCols) {
			L   =   Eigen::MatrixXd::Identity(nRows,nRows);
			R	=	get_Matrix(nRowStart, nColStart, nRows, nCols).transpose();
			// V	=	A.transpose();
			computedRank   =   nRows;
		}
		else {
			// U	=	A;
			L	=	get_Matrix(nRowStart, nColStart, nRows, nCols);
			R   =   Eigen::MatrixXd::Identity(nCols,nCols);
			computedRank   =   nCols;
		}
	}
	else {
		L   =   Eigen::MatrixXd(nRows,computedRank);
		R   =   Eigen::MatrixXd(nCols,computedRank);
		for (int j=0; j<computedRank; ++j) {
			L.col(j)    =   u[j];
			R.col(j)    =   v[j];
		}
	}
	// std::cout << "Size of row index: " << rowIndex.size() << "\n";
	// std::cout << "Size of remaining row index: " << remainingRowIndex.size() << "\n";
	return;
}

// void get_Next_Indices(int j, int k, const std::set<int>& allowed_Row_Indices, const std::set<int>& allowed_Col_Indices, int& jNew, int& kNew, double& max) {
// 	Eigen::VectorXd row	=	get_Row(j+nRowStart, nColStart, nCols);
// 	Eigen::VectorXd col	=	get_Col(k+nColStart, nRowStart, nRows);
// 	int jtemp, ktemp;
// 	Eigen::VectorXd rowtemp, coltemp;
// 	double maxRowCol, maxColRow;
// 	//	Row then column
// 	kRowCol	=	max_Abs_Vector(row, allowed_Row_Indices, maxRowCol);
// 	coltemp	=	get_Col(ktemp+nColStart, nRowStart, nRows);
// 	jRowCol	=	max_Abs_Vector(coltemp, allowed_Col_Indices, maxRowCol);
// 	//	Column then row
// 	jColRow	=	max_Abs_Vector(col, allowed_Col_Indices, maxColRow);
// 	rowtemp	=	get_Row(jtemp+nRowStart, nColStart, nCols);
// 	kColRow	=	max_Abs_Vector(rowtemp, allowed_Row_Indices, maxColRow);
// 	if (fabs(maxColRow)>fabs(maxRowCol)) {
// 		jNew	=	jColRow;
// 		kNew	=	kColRow;
// 		max		=	maxColRow;
// 	}
// 	else {
// 		jNew	=	jRowCol;
// 		kNew	=	kRowCol;
// 		max		=	maxRowCol;
// 	}
// }
//
// void HODLR_Matrix::my_Rook_Piv_HODLR_Matrix() {
// 	std::vector<int> rowIndex;		///	This stores the row indices, which have already been used.
// 	std::vector<int> colIndex;		///	This stores the column indices, which have already been used.
// 	std::set<int> remainingRowIndex;/// Remaining row indicies
// 	std::set<int> remainingColIndex;/// Remaining row indicies
// 	std::vector<Eigen::VectorXd> u;	///	Stores the column basis.
// 	std::vector<Eigen::VectorXd> v;	///	Stores the row basis.
//
// 	computedRank	=	0;
// 	for (int j=0; j<M; ++j) {
// 		remainingRowIndex.insert(j);
// 	}
// 	for (int k=0; k<N; ++k) {
// 		remainingColIndex.insert(k);
// 	}
// }

#endif /*__HODLR_Matrix__*/