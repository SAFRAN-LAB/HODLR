//
// Matrix_Entry_Routine.cpp
//
// MODIFY THIS FILE TO CHANGE THE DEFINITION
// OF YOUR SYSTEM MATRIX
//
// DO NOT MODIFY THE FUNCTION'S CALLING SEQUENCE
//
//
//  Created by Travis Askham on 5/1/2014
//
//

#include "Matrix_Entry_Routine.hpp"
#include <cmath>
#include <iostream>

using namespace std;

void Matrix_Entry_Routine(double * xyz, double * xyztarg, double * work,
			  unsigned i, unsigned j, unsigned nDim, double & val)

{
  double R2 = (xyz[nDim*j]-xyztarg[nDim*i])*(xyz[nDim*j]-xyztarg[nDim*i]);
  for (int ii = 1; ii < nDim; ii++) {
    R2 = R2 + (xyz[nDim*j+ii]-xyztarg[nDim*i+ii])*(xyz[nDim*j+ii]-xyztarg[nDim*i+ii]);
  }
  R2 = R2*work[j]*work[j];
  if ( R2 > 1.0E-15)
    val = R2*R2*log(R2)*0.5/(1+R2*R2);
  //    val = exp(-R2);
  else
    val = 0.0;
  //    val = 1.0; 
  return;
}


