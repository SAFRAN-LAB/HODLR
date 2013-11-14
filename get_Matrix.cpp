#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "get_Matrix.hpp"

using namespace std;
using namespace Eigen;


VectorXd Theta;

/****************************************************************/
/*	FUNCTION:	set_Locations				*/
/*								*/
/*	Sets the location of points in space.			*/
/****************************************************************/
void set_Locations(unsigned N) {
	Theta	=	VectorXd::Random(N);
	sort(Theta.data(), Theta.data()+Theta.size());
}


/****************************************************************/
/*	FUNCTION:   get_Matrix_Entry                         	*/
/*                                                              */
/*	Obtains an entry of the matrix                         	*/
/****************************************************************/
double get_Matrix_Entry(const unsigned i, const unsigned j) {
	double R	=	fabs(Theta(i)-Theta(j));
//		return exp(-R);
//	return sin(R)/R;
	return exp(-R*R);
	//		return exp(-sqrt(R2));
	//		return 1.0/(1.0+R2);
	
	//		double R	=	fabs(Theta(i)-Theta(j));
	//		return exp(-R);
	//		return sin(R)/R;
	//		return 1.0/R;
}


/****************************************************************/
/*	FUNCTION:   get_Matrix					*/
/*                                                              */
/*	Obtains a sub-matrix of the matrix			*/
/****************************************************************/
void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A) {
	A   =   MatrixXd(n_Rows,n_Cols);
	for (unsigned i=0; i<n_Rows; ++i) {
		for (unsigned j=0; j<n_Cols; ++j) {
			A(i,j)  =   get_Matrix_Entry(start_Row+i, start_Col+j);
		}
	}
}

/****************************************************************/
/*  FUNCTION:   get_Matrix_Row                                  */
/*                                                              */
/*  Obtains a row of a sub-matrix                               */
/****************************************************************/
void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v) {
    v   =   VectorXd(n_Cols);
    for (unsigned j=0; j<n_Cols; ++j) {
        v(j)    =   get_Matrix_Entry(row_Index, start_Col+j);
    }
}

/****************************************************************/
/*  FUNCTION:   get_Matrix_Col                                  */
/*                                                              */
/*  Obtains a column of a sub-matrix                            */
/****************************************************************/
void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v) {
    v   =   VectorXd(n_Rows);
    for (unsigned j=0; j<n_Rows; ++j) {
        v(j)    =   get_Matrix_Entry(start_Row+j, col_Index);
    }
}

/****************************************************************/
/*  FUNCTION:   max_Abs_Vector                                  */
/*                                                              */
/*  Obtains index and value of maximum of the absolute entry    */
/*  in a vector                                                 */
/****************************************************************/
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