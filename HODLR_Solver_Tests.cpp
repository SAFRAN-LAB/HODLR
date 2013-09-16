/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran
 *  \version 3.1
 */
/*! \file HODLR_Solver_Tests.cpp
 */

/*!
 *  This is a test file for the HODLR_Solver. Line 57: "HODLR_Solve(rhs, n_Leaf, fast_sol);" is the actual call to the solver; All the rest are some additional embellishments.
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
    MatrixXd fast_sol, eigen_sol, rhs;
    clock_t start, end;
    double time, error;
    unsigned N;
    
    for (int n=0; n<5; ++n) {
        N   =   n_Leaf*pow(2.0,n);
        MatrixXd exact_sol  =   MatrixXd::Random(N,1);
        
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
        
        cout << "Fast Solver..." << endl;
        start   =   clock();
        HODLR_Solve(rhs, n_Leaf, fast_sol);
        end     =   clock();
        time    =   double(end-start)/double(CLOCKS_PER_SEC);
        cout << "Time taken by the fast solver to solve the linear system is: " << time << endl << endl;
        
        cout << "Computing the relative error of the fast solver..." << endl;
        error   =   (fast_sol-exact_sol).norm()/exact_sol.norm();
        cout << "The relative error in the solution is: " << error << endl << endl;
        
        cout << "Conventional Direct Solver..." << endl;
        start       =   clock();
        eigen_sol   =   K.fullPivLu().solve(rhs);
        end         =   clock();
        time        =   double(end-start)/double(CLOCKS_PER_SEC);
        cout << "Time taken by the conventional direct solver to solve the linear system is: " << time << endl << endl;
        
        cout << "Computing the relative error of the conventional direct solver..." << endl;
        error   =   (eigen_sol-exact_sol).norm()/exact_sol.norm();
        cout << "The relative error in the solution is: " << error << endl << endl;
    }
}
