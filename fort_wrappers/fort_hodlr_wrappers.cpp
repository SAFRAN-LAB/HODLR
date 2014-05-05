//
// fort_hodlr_wrappers.cpp
//
//
//  Created by Travis Askham on 5/1/14.
//
//

#include <iostream>
#include <Eigen/Dense>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iomanip>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "HODLR_Extern.hpp"
#include "fort_hodlr_wrappers.hpp"

using namespace std;
using namespace Eigen;


void initialize_matrix_wrap(HODLR_Tree<Extern_Kernel> ** A, 
			    Extern_Kernel ** kernel, 
			    unsigned N, unsigned nDim, unsigned nLeaf, 
			    double * xyz, double * xyztarg, 
			    double * work, double * diag, double eps)
/*********************************************************************

Creates a HODLR_Tree<Extern_Kernel> at (*A) with kernel at (*kernel).


Input:

The kernel is created by calling the Extern_Kernel constructor with 
the parameters:

N - system size.
nDim - generally, the number of dimensions
xyz - real array used by Matrix_Entry_Routine
xyztarg - real array used by Matrix_Entry_Routine
work - real array used by Matrix_Entry_Routine

Other input:

diag - real array of length N containing the matrix diagonal
eps - tolerance used by assemble_Matrix 

Output:

A - the system matrix stored as a HODLR_Tree

 *********************************************************************/
{

  VectorXd diagonal(N);

  for (int i = 0; i < N; i++)
    diagonal(i) = diag[i];

  (*kernel) = new Extern_Kernel(N,nDim,xyz,xyztarg,work);
  (*A) = new HODLR_Tree<Extern_Kernel>(*kernel, N, nLeaf);

  (*A)->assemble_Matrix(diagonal, eps);

  return;
  
}

void matrix_multiply_wrap(HODLR_Tree<Extern_Kernel> ** A, double * x, 
			  double * b, unsigned N, unsigned nRow, 
			  unsigned nCol)
/*********************************************************************

Preceded by a call to initialize_matrix_wrap

calculates b = Ax

Input:

A - a HODLR_Tree corresponding to the system matrix.
x - array containing the vector (or matrix) to be multiplied
N - system size
nRow - number of rows in x (and in the returned b), i.e. if x is a 
matrix, x(i,j) in fortran is at x[(j-1)*nRow+i-1]
nCol - number of columns in x

Output:

b - array containing the resulting vector (or matrix) with same 
index convention as for x

 *********************************************************************/
{
  MatrixXd xx(N,nCol);
  MatrixXd bb(N,nCol);
  
  for (int j = 0; j < nCol; j++) {
    for (int i = 0; i < N; i ++) {
      xx(i,j) = x[nRow*j + i]; 
    }
  }
  
  (*A)->matMatProduct(xx,bb);

  for (int j = 0; j < nCol; j++) {
    for (int i = 0; i < N; i ++) {
      b[nRow*j + i] = bb(i,j); 
    }
  }

  return;

}

void matrix_factor_wrap(HODLR_Tree<Extern_Kernel> ** A)
/*********************************************************************

Preceded by a call to initialize_matrix_wrap 

factorizes the matrix A

Input:

A - a HODLR_Tree corresponding to the system matrix.

Output:

The factorized A

 *********************************************************************/

{

  (*A)->compute_Factor();

  return;

}

void matrix_solve_wrap(HODLR_Tree<Extern_Kernel> ** A, double * x, 
		       double * b, unsigned N, unsigned nRow, 
		       unsigned nCol)
/*********************************************************************

Preceded by a call to initialize_matrix_wrap and matrix_factor_wrap

calculates x = A^(-1) b

Input:

A - a HODLR_Tree corresponding to the system matrix.
b - array containing the vector (or matrix) of the right hand side
N - system size
nRow - number of rows in b (and in the returned x), i.e. if b is a 
matrix, b(i,j) in fortran is at b[(j-1)*nRow+i-1]
nCol - number of columns in b

Output:

x - array containing the resulting vector (or matrix) with same 
index convention as for b

 *********************************************************************/
{

  MatrixXd xx(N,nCol);
  MatrixXd bb(N,nCol);
  
  for (int j = 0; j < nCol; j++) {
    for (int i = 0; i < N; i ++) {
      bb(i,j) = b[nRow*j + i]; 
    }
  }

  (*A)->solve(bb,xx);

  for (int j = 0; j < nCol; j++) {
    for (int i = 0; i < N; i ++) {
      x[nRow*j + i] = xx(i,j); 
    }
  }

  return;

}

void matrix_determinant_wrap(HODLR_Tree<Extern_Kernel> ** A, 
			     double * determinant)
/*********************************************************************

Preceded by a call to initialize_matrix_wrap and matrix_factor_wrap

computes the log determinant of A

Input:

A - a HODLR_Tree corresponding to the system matrix.

Output:

deteriminant - the log determinant of the matrix

 *********************************************************************/

{

  (*A)->compute_Determinant(*determinant);

}
