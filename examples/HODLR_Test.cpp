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

#include "HODLR_Tree.hpp"
#include "HODLR_Kernel.hpp"
#include "KDTree.hpp"

using namespace std;
using namespace Eigen;

class Test_Kernel : public HODLR_Kernel {

private:

#ifdef ONE
	VectorXd Theta;
	const static unsigned nDim	=	1;
#elif TWO
	MatrixXd Theta;
	const static unsigned nDim	=	2;
#elif THREE
	MatrixXd Theta;
	const static unsigned nDim	=	3;
#endif

public:
	Test_Kernel (unsigned N) {
#ifdef ONE
		Theta	=        VectorXd::Random(N);
		sort(Theta.data(), Theta.data()+Theta.size());
#else
		Theta	=	MatrixXd::Random(N, nDim);
		get_KDTree_Sorted(Theta,0);
#endif
	};

	double get_Matrix_Entry(const unsigned i, const unsigned j) {
#ifdef ONE
		double R	=	fabs(Theta(i)-Theta(j));
		#ifdef	GAUSSIAN
			return exp(-R*R);
		#elif	EXPONENTIAL
			return exp(-R);
		#elif	SINC
			return sin(R)/R;
		#elif	QUADRIC
			return 1.0+R*R;
		#elif	INVERSEQUADRIC
			return 1.0/(1.0+R*R);
		#elif	MULTIQUADRIC
			return sqrt(1.0+R*R);
		#elif	INVERSEMULTIQUADRIC
			return 1.0/sqrt(1.0+R*R);
		#elif	R2LOGR
				return R*R*log(R);
		#elif	LOGR
				return log(R);
		#elif	ONEOVERR
				return 1.0/R;
		#elif	LOG1R
				return log(1+R);
		#endif
#else
		double R2	=	(Theta(i,0)-Theta(j,0))*(Theta(i,0)-Theta(j,0));
		for (unsigned k=1; k<nDim; ++k) {
			R2	=	R2+(Theta(i,k)-Theta(j,k))*(Theta(i,k)-Theta(j,k));
		}
		#ifdef	GAUSSIAN
			return exp(-R2);
		#elif	EXPONENTIAL
			return exp(-sqrt(R2));
		#elif	SINC
			double R	=	sqrt(R2);
			return sin(R)/R;
		#elif	QUADRIC
			return 1.0+R2;
		#elif	INVERSEQUADRIC
			return 1.0/(1.0+R2);
		#elif	MULTIQUADRIC
			return sqrt(1.0+R2);
		#elif	INVERSEMULTIQUADRIC
			return 1.0/sqrt(1.0+R2);
		#elif	R2LOGR
				return 0.5*R2*log(R2);
		#elif	LOGR
				return 0.5*log(R2);
		#elif	ONEOVERR
				return 1.0/sqrt(R2);
		#elif	LOG1R
				return log(1+sqrt(R2));
		#endif
#endif
	};
};

int main() {
	srand (time(NULL));

	unsigned N	=	100000;
	unsigned nRhs	=	1;
	unsigned nLeaf	=	50;
	double tolerance=	1e-14;

	Test_Kernel kernel(N);

	MatrixXd x	=	MatrixXd::Random(N, nRhs);
	MatrixXd b, xSol;

	cout << endl << "Number of particles is: " << N << endl;
	clock_t start, end;

	cout << endl << "Setting things up..." << endl;
	start	=	clock();
	HODLR_Tree<Test_Kernel>* A	=	new HODLR_Tree<Test_Kernel>(&kernel, N, nLeaf);
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
