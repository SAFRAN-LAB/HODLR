#include <Eigen/Dense>
#include <vector>
#include <cstdlib>
#include "partial_Piv_LU.hpp"
#include "get_Matrix.hpp"

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

void partial_Piv_LU(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const double tolerance, unsigned& computed_Rank, MatrixXd& U, MatrixXd& V) {
    
	//  If the matrix is small enough, do not do anything
	unsigned tolerable_Rank =   5;
    if (n_Cols <= tolerable_Rank){
        get_Matrix(start_Row, start_Col, n_Rows, n_Cols, U);
        V               =   MatrixXd::Identity(n_Cols, n_Cols);
        computed_Rank   =   n_Cols;
        return;
	}
    else if (n_Rows <= tolerable_Rank){
        U               =   MatrixXd::Identity(n_Rows, n_Rows);
        get_Matrix(start_Row, start_Col, n_Rows, n_Cols, V);
        computed_Rank   =   n_Rows;
        return;
    }

    vector<int> rowIndex;	//	This stores the row indices, which have already been used.
    vector<int> colIndex;	//	This stores the column indices, which have already been used.
    vector<VectorXd> u;		//	Stores the column basis.
    vector<VectorXd> v;		//	Stores the row basis.

	srand (time(NULL));
    double max, Gamma, unused_max;
    
    /*  INITIALIZATION  */
    
    //  Initialize the matrix norm and the the first row index
    double matrix_Norm  =   0;
    rowIndex.push_back(0);
    
    unsigned pivot;
    
    computed_Rank   =   0;
    
    VectorXd a, row, col;
    
    double row_Squared_Norm, row_Norm, col_Squared_Norm, col_Norm;

    //  Repeat till the desired tolerance is obtained
    do {
        //  Generation of the row
        get_Matrix_Row(start_Col, n_Cols, start_Row+rowIndex.back(), a);
        //  Row of the residuum and the pivot column
        row =   a;
        for (unsigned l=0; l<computed_Rank; ++l) {
            row =   row-u[l](rowIndex.back())*v[l];
        }

        pivot   =   max_Abs_Vector(row, colIndex, max);
        
        unsigned max_tries  =   50;
        unsigned count      =   0;
        unsigned count1     =   0;
        
        //  This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows.
        while (fabs(max)<tolerance && count < max_tries) {
            unsigned new_rowIndex;
            rowIndex.pop_back();
            do {
                new_rowIndex   =   rand()%n_Rows;
                ++count1;
            } while (find(rowIndex.begin(),rowIndex.end(),new_rowIndex)!=rowIndex.end() && count1 < max_tries);
            count1  =   0;
            rowIndex.push_back(new_rowIndex);
            
            //  Generation of the row
            get_Matrix_Row(start_Col, n_Cols, start_Row+rowIndex.back(), a);
            
            //  Row of the residuum and the pivot column
            row =   a;
            for (unsigned l=0; l<computed_Rank; ++l) {
                row =   row-u[l](rowIndex.back())*v[l];
            }
            pivot   =   max_Abs_Vector(row, colIndex, max);
            ++count;
        }

        if (count == max_tries) break;
        
        count = 0;
        
        colIndex.push_back(pivot);
        
        //  Normalizing constant
        Gamma   =   1.0/max;
        
        //  Generation of the column
        get_Matrix_Col(start_Row, n_Rows, start_Col+colIndex.back(), a);
        
        //  Column of the residuum and the pivot row
        col =   a;
        for (unsigned l=0; l<computed_Rank; ++l) {
            col =   col-v[l](colIndex.back())*u[l];
        }
        pivot   =   max_Abs_Vector(col, rowIndex, unused_max);
        
        //  This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
        while (fabs(max)<tolerance && count < max_tries) {
            colIndex.pop_back();
            unsigned new_colIndex;
            do {
                new_colIndex   =   rand()%n_Cols;
            } while (find(colIndex.begin(),colIndex.end(),new_colIndex)!=colIndex.end() && count1 < max_tries);
            count1  =   0;
            colIndex.push_back(new_colIndex);
            
            //  Generation of the column
            get_Matrix_Col(start_Row, n_Rows, start_Col+colIndex.back(), a);
            
            //  Column of the residuum and the pivot row
            col =   a;
            for (unsigned l=0; l<computed_Rank; ++l) {
                col =   col-u[l](colIndex.back())*v[l];
            }
            pivot   =   max_Abs_Vector(col, rowIndex, unused_max);
            ++count;
        }
        
        if (count == max_tries) break;
        
        count = 0;
        
        rowIndex.push_back(pivot);
        
        //  New vectors
        u.push_back(Gamma*col);
        v.push_back(row);
        
        //  New approximation of matrix norm
        row_Squared_Norm    =   row.squaredNorm();
        row_Norm            =   sqrt(row_Squared_Norm);
        
        col_Squared_Norm    =   col.squaredNorm();
        col_Norm            =   sqrt(col_Squared_Norm);
        
        matrix_Norm         =   matrix_Norm +   Gamma*Gamma*row_Squared_Norm*col_Squared_Norm;
        
        for (unsigned j=0; j<computed_Rank; ++j) {
            matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
        }
        ++computed_Rank;
    } while (row_Norm*col_Norm > fabs(max)*tolerance*matrix_Norm && computed_Rank <= fmin(n_Rows, n_Cols));

    //  If the computed_Rank is close to full-rank then return the trivial full-rank decomposition
    if (computed_Rank>=fmin(n_Rows, n_Cols)) {
        if (n_Rows < n_Cols) {
            U   =   MatrixXd::Identity(n_Rows,n_Rows);
            get_Matrix(start_Row, start_Col, n_Rows, n_Cols, V);
            computed_Rank   =   n_Rows;
            return;
        }
        else {
            get_Matrix(start_Row, start_Col, n_Rows, n_Cols, U);
            V   =   MatrixXd::Identity(n_Cols,n_Cols);
            computed_Rank   =   n_Cols;
            return;
        }
    }

    U   =   MatrixXd(n_Rows,computed_Rank);
    V   =   MatrixXd(computed_Rank,n_Cols);
    for (unsigned j=0; j<computed_Rank; ++j) {
        U.col(j)    =   u[j];
        V.row(j)    =   v[j];
    }
}