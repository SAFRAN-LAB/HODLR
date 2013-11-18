#ifndef __get_Matrix_hpp__
#define __get_Matrix_hpp__


#include<vector>
#include<cmath>
using namespace std;

#include"Eigen/Dense"
using namespace Eigen;


/****************************************************************/
/*	FUNCTION:	set_Locations				*/
/*								*/
/*	Sets the location of points in space.			*/
/****************************************************************/
void set_Locations(unsigned N);

/****************************************************************/
/*	FUNCTION:   get_Matrix_Entry                         	*/
/*                                                              */
/*	Obtains an entry of the matrix                         	*/
/****************************************************************/
double get_Matrix_Entry(const unsigned i, const unsigned j);


/****************************************************************/
/*	FUNCTION:   get_Matrix					*/
/*                                                              */
/*	Obtains a sub-matrix of the matrix			*/
/****************************************************************/
void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A);

/****************************************************************/
/*  FUNCTION:   get_Matrix_Row                                  */
/*                                                              */
/*  Obtains a row of a sub-matrix                               */
/****************************************************************/
void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v);

/****************************************************************/
/*  FUNCTION:   get_Matrix_Col                                  */
/*                                                              */
/*  Obtains a column of a sub-matrix                            */
/****************************************************************/
void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v);

/****************************************************************/
/*  FUNCTION:   max_Abs_Vector                                  */
/*                                                              */
/*  Obtains index and value of maximum of the absolute entry    */
/*  in a vector                                                 */
/****************************************************************/
unsigned max_Abs_Vector(const VectorXd& v, const vector<int>& not_Allowed_Indices, double& max);

#endif