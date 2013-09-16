/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran
 *  \version 3.1
 */
/*! \file HODLR_Solver.hpp
 */

/*!
*****************************************
*   PREREQUISITES FOR USING THIS HEADER *
*****************************************

The main and only prerequisite to use this header is to have the linear algebra package EIGEN available at: http://eigen.tuxfamily.org/index.php?title=Main_Page.
*/

#include "iostream"
#include "Eigen/Dense"
#include "vector"
#include "cmath"
#include "ctime"

using namespace std;
using namespace Eigen;

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      HODLR_Solve(const MatrixXd& b, const unsigned n_Leaf, MatrixXd& x)  */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      User calls this function to perform the fast solve                  */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      b       -   Right hand side;            Type:   MatrixXd            */
/*      n_Leaf  -   Size of the mininimum sub-matrix that can be solved     */
/*                  without the fast solver;    Type:   unsigned            */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      x       -   Solution;           Type:   MatrixXd                    */
/*                                                                          */
/****************************************************************************/
void HODLR_Solve(const MatrixXd& b, const unsigned n_Leaf, MatrixXd& x);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      fast_Solve(const unsigned start_Row, const unsigned start_Col,      */
/*          const unsigned n_Rows, const unsigned n_Cols,                   */
/*          const unsigned n_Leaf, const MatrixXd& b, MatrixXd& x)          */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Calls itself recursively to perform fast solve on smaller matrices  */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      start_Row   -   Starting row of the sub-matrix;     Type: unsigned  */
/*      start_Col   -   Starting column of the sub-matrix;  Type: unsigned  */
/*      n_Rows      -   Number of rows of the sub-matrix;   Type: unsigned  */
/*      n_Cols      -   Number of columns of the sub-matrix;Type: unsigned  */
/*      n_Leaf      -   Size of the mininimum sub-matrix that can be solved */
/*                          without the fast solver;        Type: unsigned  */
/*      b           -   Right hand side;                    Type: MatrixXd  */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      x           -   Solution;                           Type: MatrixXd  */
/*                                                                          */
/****************************************************************************/
void fast_Solve(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const unsigned n_Leaf, const MatrixXd& b, MatrixXd& x);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      partial_Piv_LU(const unsigned start_Row, const unsigned start_Col,  */
/*          const unsigned n_Rows, const unsigned n_Cols,                   */
/*          const double tolerance, unsigned& computed_Rank, MatrixXd& U,   */
/*          MatrixXd& V)                                                    */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains low-rank decomposition to desired tolerance of a sub-matrix */
/*          using partial pivoted LU algorithm.                             */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      start_Row       -   Starting row of the sub-matrix;                 */
/*                                                          Type: unsigned  */
/*                                                                          */
/*      start_Col       -   Starting column of the sub-matrix;              */
/*                                                          Type: unsigned  */
/*                                                                          */
/*      n_Rows          -   Number of rows of the sub-matrix;               */
/*                                                          Type: unsigned  */
/*                                                                          */
/*      n_Cols          -   Number of columns of the sub-matrix;            */
/*                                                          Type: unsigned  */
/*                                                                          */
/*      tolerance       -   Tolerance of low-rank approximation;            */
/*                                                          Type: double    */
/*                                                                          */
/*      computed_Rank   -   Rank of the decomposition;                      */
/*                                                          Type: unsigned  */
/*                                                                          */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      Obtains the low-rank decomposition as UV                            */
/*      U   -   Column basis                                                */
/*      V   -   Row basis                                                   */
/*                                                                          */
/****************************************************************************/
void partial_Piv_LU(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const double tolerance, unsigned& computed_Rank, MatrixXd& U, MatrixXd& V);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      assemble_Matrix(const unsigned M, const unsigned N)                 */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Assembles the matrix of size M by N                                 */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      M   -   Number of rows;     Type:   unsigned                        */
/*      N   -   Number of columns;  Type:   unsigned                        */
/*                                                                          */
/****************************************************************************/
void assemble_Matrix(const unsigned M, const unsigned N);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      double get_Matrix_Entry(const unsigned i, const unsigned j)         */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains an entry of the matrix                                      */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      i       -   Row of the matrix;      Type:   unsigned                */
/*      j       -   Column of the matrix;   Type:   unsigned                */
/*                                                                          */
/*  RETURNS:                                                                */
/*      Returns the (i,j)^th entry of the matrix of type double             */
/*                                                                          */
/****************************************************************************/
double get_Matrix_Entry(const unsigned i, const unsigned j);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      get_Matrix(const unsigned start_Row, const unsigned start_Col,      */
/*          const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A)      */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains a sub-matrix of the matrix                                  */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      start_Row   -   Starting row of the sub-matrix;     Type: unsigned  */
/*      start_Col   -   Starting column of the sub-matrix;  Type: unsigned  */
/*      n_Rows      -   Number of rows of the sub-matrix;   Type: unsigned  */
/*      n_Cols      -   Number of columns of the sub-matrix;Type: unsigned  */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      A           -   The desired sub-matrix              Type: MatrixXd  */
/*                                                                          */
/****************************************************************************/
void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols,     */
/*              const unsigned row_Index, VectorXd& v)                      */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains a row of the matrix                                         */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      row_Index   -   Desired row of the matrix;          Type: unsigned  */
/*      start_Col   -   Starting column of the desired row; Type: unsigned  */
/*      n_Cols      -   Number of columns of desired row;   Type: unsigned  */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      v           -   The desired row                     Type: VectorXd  */
/*                                                                          */
/****************************************************************************/
void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows,     */
/*              const unsigned col_Index, VectorXd& v)                      */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains a column of the matrix                                      */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      col_Index   -   Desired column of the matrix;       Type: unsigned  */
/*      start_Row   -   Starting row of the desired column; Type: unsigned  */
/*      n_Rows      -   Number of rows of desired column;   Type: unsigned  */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      v           -   The desired column                  Type: VectorXd  */
/*                                                                          */
/****************************************************************************/
void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v);

/****************************************************************************/
/*                                                                          */
/*  FUNCTION:                                                               */
/*      unsigned max_Abs_Vector(const VectorXd& v, double& max,             */
/*              const vector<int>& not_Allowed_Indices)                     */
/*                                                                          */
/*  PURPOSE OF THIS FUNCTION:                                               */
/*      Obtains index & value of the maximum absolute entry in the vector   */
/*                                                                          */
/*  INPUTS:                                                                 */
/*      v   -   Vector whose maximum absolute entry is desired;             */
/*                                                      Type: unsigned      */
/*                                                                          */
/*      not_Allowed_Indices -   Indices which are not allowed;              */
/*                                                      Type: vector<int>   */
/*                                                                          */
/*  OUTPUTS:                                                                */
/*      max -   The maximum absolute entry in the vector                    */
/*                                                                          */
/*  RETURNS:                                                                */
/*      Returns index of the maximum absolute entry in the vector           */
/*                                                                          */
/****************************************************************************/
unsigned max_Abs_Vector(const VectorXd& v, double& max, const vector<int>& not_Allowed_Indices);