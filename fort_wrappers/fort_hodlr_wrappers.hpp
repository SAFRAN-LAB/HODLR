//
// fort_hodlr_wrappers.hpp
//
//
//  Created by Travis Askham on 5/1/14.
//
//
#ifndef __FORT_HODLR_WRAPPERS__
#define __FORT_HODLR_WRAPPERS__

#include "HODLR_Tree.hpp"
#include "HODLR_Extern.hpp"

extern "C"
{

  void initialize_matrix_wrap(HODLR_Tree<Extern_Kernel> ** A, 
			      Extern_Kernel ** kernel, 
			      unsigned N, unsigned nDim, unsigned nLeaf, 
			      double * xyz, double * xyztarg, double * work, 
			      double * diag, double eps);

  void matrix_multiply_wrap(HODLR_Tree<Extern_Kernel> ** A, 
			    double * x, double * b, unsigned N, 
			    unsigned nRow, unsigned nCol);
  
  void matrix_factor_wrap(HODLR_Tree<Extern_Kernel> ** A);

  void matrix_solve_wrap(HODLR_Tree<Extern_Kernel> ** A, double * x, 
			 double * b, unsigned N, unsigned nRow, 
			 unsigned nCol);

  void matrix_determinant_wrap(HODLR_Tree<Extern_Kernel> ** A, 
			       double * determinant);
  
}


#endif /* (defined(__FORT_HODLR_WRAPPERS__) */
