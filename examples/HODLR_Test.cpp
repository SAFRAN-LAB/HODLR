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
#include <vector>
#include <iomanip>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

using std::cout;
using std::cout;
using std::vector;
using std::setprecision;
using std::sort;
using namespace Eigen;

class Test_Kernel : public HODLR_Matrix {

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

	unsigned N	=	50000;
	unsigned nRhs	=	1;
	unsigned nLeaf	=	100;
	double tolerance=	1e-15;

	Test_Kernel kernel(N);

	MatrixXd xExact	=	MatrixXd::Random(N, nRhs);
	MatrixXd bExact(N,nRhs), bFast(N,nRhs), xFast(N,nRhs);

	cout << endl << "Number of particles is: " << N << endl;
	clock_t start, end;

	cout << endl << "Setting things up..." << endl;
	start	=	clock();
	HODLR_Tree<Test_Kernel>* A	=	new HODLR_Tree<Test_Kernel>(&kernel, N, nLeaf);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Assembling the matrix in HODLR form..." << endl;
	start			=	clock();
	VectorXd diagonal	=	2.0*VectorXd::Ones(N)+VectorXd::Random(N);
	A->assemble_Matrix(diagonal, tolerance, 's');
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

        cout << endl << "Exact matrix vector product..." << endl;
        start           =       clock();
        for (int i=0; i<N; ++i) {
                bExact(i,0)             =       diagonal(i)*xExact(i,0);
                for (int j=0; j<i; ++j) {
                        bExact(i,0)     =       bExact(i,0)+kernel.get_Matrix_Entry(i, j)*xExact(j,0);
                }
                for (int j=i+1; j<N; ++j) {
                        bExact(i,0)     =       bExact(i,0)+kernel.get_Matrix_Entry(i, j)*xExact(j,0);
                }
        }
        end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Fast matrix matrix product..." << endl;
	start		=	clock();
	A->matMatProduct(xExact, bFast);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Factoring the matrix..." << endl;
	start		=	clock();
	A->compute_Factor();
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Solving the system..." << endl;
	start		=	clock();
	A->solve(bExact, xFast);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Error in computed solution: " << (xFast-xExact).norm()/xExact.norm()<< endl;

	cout << endl << "Error in matrix matrix product: " << (bFast-bExact).cwiseAbs().maxCoeff() << endl;
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
