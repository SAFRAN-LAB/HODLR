/*!
 \class HODLR_Matrix

 \brief This class is for HODLR matrix.

 \note

 \author $Dan Foreman-Mackey$

 \version

 \date $Date: January 31st, 2014$
 */

#ifndef __HODLR_MATRIX_HPP__
#define __HODLR_MATRIX_HPP__

#include <vector>
#include <cmath>
#include <Eigen/Dense>

using std::vector;
using namespace Eigen;

class HODLR_Matrix {

public:
	HODLR_Matrix () {};

        /*!
         Allows access to the matrix elements, i.e., returns the (i,j)th element of the HODLR matrix.
         */
	virtual double get_Matrix_Entry(const unsigned i, const unsigned j) {
		return 0.0;
	};

        /*!
         Allows access to the sub-matrix, i.e., returns a 'n_Rows' by 'n_Cols' sub-matrix starting at the index (start_Row, start_Column) of the HODLR matrix.
         */
	void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A) {
		A   =   MatrixXd(n_Rows,n_Cols);
		for (unsigned i=0; i<n_Rows; ++i) {
			for (unsigned j=0; j<n_Cols; ++j) {
				A(i,j)  =   get_Matrix_Entry(start_Row+i, start_Col+j);
			}
		}
	};

        /*!
         Allows access to the rows of the HODLR matrix, i.e., returns the 'row_Index'th row with the column starting at 'start_Col' and ending at 'start_Col+n_Cols'.
         */
	void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v) {
		v   =   VectorXd(n_Cols);
		for (unsigned j=0; j<n_Cols; ++j) {
			v(j)    =   get_Matrix_Entry(row_Index, start_Col+j);
		}
	}

        /*!
         Allows access to the columns of the HODLR matrix, i.e., returns the 'col_Index'th column with the row starting at 'start_Row' and ending at 'start_Row+n_Rows'.
         */
	void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v) {
		v   =   VectorXd(n_Rows);
		for (unsigned j=0; j<n_Rows; ++j) {
			v(j)    =   get_Matrix_Entry(start_Row+j, col_Index);
		}
	}

        /*!
         Returns the maximum absolute value and the index of maximum element of the vector.
         */
	unsigned max_Abs_Vector(const VectorXd& v, const vector<int>& not_Allowed_Indices, double& max) {
		unsigned n      =   v.size();
		max             =   v(0);
		unsigned index  =   0;
		for(unsigned j=0; j<n; ++j){
			if(find(not_Allowed_Indices.begin(),not_Allowed_Indices.end(),j)==not_Allowed_Indices.end()) {
				if(fabs(v(j))>fabs(max)){
					max     =   v(j);
					index   =   j;
				}
			}
		}
		return index;
	}

};

#endif /* defined(__HODLR_MATRIX_HPP__) */
