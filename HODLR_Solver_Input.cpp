/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran
 *  \version 3.1
 */
/*! \file HODLR_Solver_Input.cpp
 */

/*!
 *  This is a sample input file for the HODLR_Solver. Line 55: "HODLR_Solve(rhs, n_Leaf, fast_sol);" is the actual call to the solver; All the rest are some additional embellishments.
 */

#include "iostream"
#include "Eigen/Dense"
#include "vector"
#include "cmath"
#include "ctime"
#include "HODLR_Solver.hpp"

using namespace std;
using namespace Eigen;

extern MatrixXd K;

int main(){
    srand (time(NULL));
    unsigned n_Leaf     =   100;
    unsigned N          =   8000;
    MatrixXd exact_sol  =   MatrixXd::Random(N,1);
    
    MatrixXd fast_sol, rhs;
    clock_t start, end;
    double time, error;
    
    cout << "Size of the matrix is: " << N << endl << endl;
    
    cout << "Assembling the matrix..." << endl;
    start   =   clock();
    assemble_Matrix(N, N);
    end     =   clock();
    time    =   double(end-start)/double(CLOCKS_PER_SEC);
    cout << "Time taken to assemble the matrix is: " << time << endl << endl;

    cout << "Obtaining the rhs..." << endl;
    start   =   clock();
    rhs     =   K*exact_sol;
    end     =   clock();
    time    =   double(end-start)/double(CLOCKS_PER_SEC);
    cout << "Time taken to obtain the rhs is: " << time << endl << endl;
    
    cout << "Solving the linear system..." << endl;
    start   =   clock();
    HODLR_Solve(rhs, n_Leaf, fast_sol);
    end     =   clock();
    time    =   double(end-start)/double(CLOCKS_PER_SEC);
    cout << "Time taken to solve the linear system is: " << time << endl << endl;
    
    cout << "Computing the relative error..." << endl;
    error   =   (fast_sol-exact_sol).norm()/exact_sol.norm();
    cout << "The relative error in the solution is: " << error << endl << endl;
}
