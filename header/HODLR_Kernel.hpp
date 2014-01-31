//
//  HODLR_Kernel.hpp
//
//
//  Created by Dan Foreman-Mackey on 2014-01-31.
//
//

#ifndef __HODLR_KERNEL_HPP__
#define __HODLR_KERNEL_HPP__

#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "KDTree.hpp"

using namespace std;
using namespace Eigen;

class HODLR_Kernel {

public:
	HODLR_Kernel () {};

	virtual double get_Matrix_Entry(const unsigned i, const unsigned j) {
		return 0.0;
	};

	void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A) {
		A   =   MatrixXd(n_Rows,n_Cols);
		for (unsigned i=0; i<n_Rows; ++i) {
			for (unsigned j=0; j<n_Cols; ++j) {
				A(i,j)  =   get_Matrix_Entry(start_Row+i, start_Col+j);
			}
		}
	};

	void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v) {
		v   =   VectorXd(n_Cols);
		for (unsigned j=0; j<n_Cols; ++j) {
			v(j)    =   get_Matrix_Entry(row_Index, start_Col+j);
		}
	}

	void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v) {
		v   =   VectorXd(n_Rows);
		for (unsigned j=0; j<n_Rows; ++j) {
			v(j)    =   get_Matrix_Entry(start_Row+j, col_Index);
		}
	}

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

#endif /* defined(__HODLR_KERNEL_HPP__) */
