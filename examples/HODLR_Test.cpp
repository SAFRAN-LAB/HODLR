//
//  HODLR_Test.cpp
//  
//
//  Created by Sivaram Ambikasaran on 11/8/13.
//
//

#include <iostream>
#include <Eigen/Dense>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include "get_Matrix.hpp"
#include "HODLR_Tree.hpp"


using namespace std;
using namespace Eigen;

int main() {
    srand (time(NULL));
	
	unsigned N	=	100000;
	unsigned nRhs	=	1;
	unsigned nLeaf	=	50;
	double tolerance=	1e-14;

	MatrixXd x	=	MatrixXd::Random(N, nRhs);
	MatrixXd b, xSol;

	cout << endl << "Number of particles is: " << N << endl;
	clock_t start, end;

	cout << endl << "Setting things up..." << endl;
	start		=	clock();
	HODLR_Tree* A	=	new HODLR_Tree(N, nLeaf);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Assembling the matrix in HODLR form..." << endl;
	start			=	clock();
	VectorXd diagonal	=	2.0*VectorXd::Ones(N);
	A->assemble_Matrix(diagonal, tolerance);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Matrix matrix product..." << endl;
	start		=	clock();
	A->matMatProduct(x, b);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Factoring the matrix..." << endl;
	start		=	clock();
	A->compute_Factor();
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Solving the system..." << endl;
	start		=	clock();
	A->solve(b, xSol);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Error in computed solution: " << (xSol-x).cwiseAbs().maxCoeff() << endl;

//	MatrixXd B;
//	cout << endl << "Assembling the entire matrix..." << endl;
//	start			=	clock();
//	get_Matrix(0, 0, N, N, B);
//	end				=	clock();
//	cout << endl << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;
//
//	cout << endl << "Exact determinant is: " << setprecision(16) << log(fabs(B.partialPivLu().determinant())) << endl;

	double determinant;
	cout << endl << "Computing the log determinant..." << endl;
	start		=	clock();
	A->compute_Determinant(determinant);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Log determinant is: " << setprecision(16) << determinant << endl;
	//
	// cout << endl << "Exact matrix matrix product..." << endl;
	// start			=	clock();
	// MatrixXd bExact	=	B*x;
	// end				=	clock();
	// cout << endl << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;
	// 
	// cout << endl << (bExact-b).cwiseAbs().maxCoeff() << endl;
}
