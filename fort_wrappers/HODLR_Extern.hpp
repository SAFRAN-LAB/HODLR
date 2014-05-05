//
//  HODLR_Extern.hpp
//
//  Created by Travis Askham on 4/30/2014
//
//  Based on HODLR_Test.cpp by Sivaram Ambikasaran
//
//  Extern_Kernel is a classed derived from HODLR_Matrix. 
//  Extern_Kernel allows construction of a kernel
//  whose matrix entries are calculated by the routine
//  Matrix_Entry_Routine (separate file) which is passed
//  various arrays stored in the Extern_Kernel class
//
//
#ifndef __HODLR_EXTERN_HPP__
#define __HODLR_EXTERN_HPP__

#include <Eigen/Dense>
#include <cstdlib>
#include <vector>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "Matrix_Entry_Routine.hpp"

using std::vector;
using std::setprecision;
using std::sort;
using namespace Eigen;

class Extern_Kernel : public HODLR_Matrix {

private:

  // some arrays used to compute matrix entries
  double * xyz;
  double * xyztarg;
  double * work;
  // dimension of ambient space
  unsigned nDim;
  // system size
  unsigned N;

public:

  Extern_Kernel (unsigned Ntemp, unsigned nDimtemp, double * xyztemp, 
		 double * xyztargtemp, double * worktemp) {
    // constructor simply assigns all of the private variables
    xyz = xyztemp;
    xyztarg = xyztargtemp;
    work = worktemp;
    nDim = nDimtemp;
    N = Ntemp;

    return;

  };

  double get_Matrix_Entry(const unsigned i, const unsigned j) {
    // to get a matrix entry, call the external routine 
    // Matrix_Entry_Routine (defined in Matrix_Entry_Routine.cpp/.hpp)
    double val;

    Matrix_Entry_Routine(xyz,xyztarg,work,i,j,nDim,val);
    
    return val;
    
  };



};



#endif /* (defined(__HODLR_EXTERN_HPP__) */
