#ifndef __partial_Piv_LU_hpp__
#define __partial_Piv_LU_hpp__

# include "Eigen/Dense"

using namespace Eigen;
using namespace std;

/********************************/
/*	PURPOSE OF EXISTENCE	*/
/********************************/

//	Obtains the low-rank decomposition of the matrix to a desired tolerance using the partial pivoting LU algorithm, i.e., given a sub-matrix 'A' and tolerance 'epsilon', computes matrices 'U' and 'V' such that ||A-UV||_F < epsilon.

/************************/
/*	INPUTS          */
/************************/

//	start_Row	-	Starting row of the sub-matrix.

//	start_Col	-	Starting column of the sub-matrix.

//	n_Rows		-	Number of rows of the sub-matrix.

//	n_Cols		-	Number of columns of the sub-matrix.

//	tolerance	-	Tolerance of low-rank approximation.

/************************/
/*	OUTPUTS		*/
/************************/

//	computed_Rank	-	Rank obtained for the given tolerance.

//	U		-	Matrix forming the column basis.

//	V		-	Matrix forming the row basis.

void partial_Piv_LU(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const double tolerance, unsigned& computed_Rank, MatrixXd& U, MatrixXd& V);

#endif