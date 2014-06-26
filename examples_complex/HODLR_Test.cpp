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
#include <complex>

#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

using std::cout;
using std::cout;
using std::vector;
using std::setprecision;
using std::sort;
using std::complex;
using namespace Eigen;

complex<double> compute_Determinant(MatrixXcd& K) {
        FullPivLU<MatrixXcd> Kinverse;
        Kinverse.compute(K);
        complex<double> determinant;
        if (K.rows()>0) {        //      Check needed when the matrix is predominantly diagonal.
                MatrixXcd LU    =       Kinverse.matrixLU();
                determinant     =       log(LU(0,0));
                for (int k=1; k<K.rows(); ++k) {
                        determinant+=log(LU(k,k));
                }
                //              Previous version which had some underflow.
                //              determinant	=	log(abs(K.determinant()));
        }
        return determinant;
};

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

	std::complex<double> get_Matrix_Entry(const unsigned i, const unsigned j) {
#ifdef ONE
		double R	=	fabs(Theta(i)-Theta(j));
		#ifdef	GAUSSIAN
                        return complex<double> (exp(-R*R),1.0);
		#elif	EXPONENTIAL
			return complex<double> (exp(-R),1.0);
		#elif	SINC
			return complex<double> (sin(R)/R,1.0);
		#elif	QUADRIC
			return complex<double> (1.0+R*R,1.0);
		#elif	INVERSEQUADRIC
			return complex<double> (1.0/(1.0+R*R),1.0);
		#elif	MULTIQUADRIC
			return complex<double> (sqrt(1.0+R*R),1.0);
		#elif	INVERSEMULTIQUADRIC
			return complex<double> (1.0/sqrt(1.0+R*R),1.0);
		#elif	R2LOGR
                        return complex<double> (R*R*log(R),1.0);
		#elif	LOGR
			return complex<double> (log(R),1.0);
		#elif	ONEOVERR
			return complex<double> (1.0/R,1.0);
		#elif	LOG1R
			return complex<double> (log(1+R),1.0);
		#endif
#else
		double R2	=	(Theta(i,0)-Theta(j,0))*(Theta(i,0)-Theta(j,0));
		for (unsigned k=1; k<nDim; ++k) {
			R2	=	R2+(Theta(i,k)-Theta(j,k))*(Theta(i,k)-Theta(j,k));
		}
		#ifdef	GAUSSIAN
			return complex<double> (exp(-R2),2.0);
		#elif	EXPONENTIAL
			return complex<double> (exp(-sqrt(R2)),2.0);
		#elif	SINC
			double R	=	sqrt(R2);
			return complex<double> (sin(R)/R,2.0);
		#elif	QUADRIC
			return complex<double> (1.0+R2,2.0);
		#elif	INVERSEQUADRIC
			return complex<double> (1.0/(1.0+R2),2.0);
		#elif	MULTIQUADRIC
			return complex<double> (sqrt(1.0+R2),2.0);
		#elif	INVERSEMULTIQUADRIC
			return complex<double> (1.0/sqrt(1.0+R2),2.0);
		#elif	R2LOGR
                        return complex<double> (0.5*R2*log(R2),2.0);
		#elif	LOGR
			return complex<double> (0.5*log(R2),2.0);
		#elif	ONEOVERR
			return complex<double> (1.0/sqrt(R2),2.0);
		#elif	LOG1R
			return complex<double> (log(1+sqrt(R2)),2.0);
		#endif
#endif
	};
};

int main() {
	srand (time(NULL));

	unsigned N	=	2000;
	unsigned nRhs	=	1;
	unsigned nLeaf	=	100;
	double tolerance=	1e-15;

	Test_Kernel kernel(N);

	MatrixXcd xExact	=	MatrixXcd::Random(N, nRhs);
	MatrixXcd bExact(N,nRhs), bFast(N,nRhs), xFast(N,nRhs);

	cout << endl << "Number of particles is: " << N << endl;
	clock_t start, end;

	cout << endl << "Setting things up..." << endl;
	start	=	clock();
	HODLR_Tree<Test_Kernel>* A	=	new HODLR_Tree<Test_Kernel>(&kernel, N, nLeaf);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Assembling the matrix in HODLR form..." << endl;
	start			=	clock();
	VectorXcd diagonal	=	4.0*VectorXcd::Ones(N);
	A->assemble_Matrix(diagonal, tolerance);
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
	//	MatrixXcd B;
	//	cout << endl << "Assembling the entire matrix..." << endl;
	//	start			=	clock();
	//	get_Matrix(0, 0, N, N, B);
	//	end				=	clock();
	//	cout << endl << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;
	//
	//	cout << endl << "Exact determinant is: " << setprecision(16) << log(fabs(B.partialPivLu().determinant())) << endl;

	complex<double> determinant;
	cout << endl << "Computing the log determinant..." << endl;
	start		=	clock();
	A->compute_Determinant(determinant);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

	cout << endl << "Log determinant is: " << setprecision(16) << determinant << endl;


        MatrixXcd K;
        kernel.get_Matrix(0, 0, N, N, K);
        for (int k=0; k<N; ++k) {
                K(k,k)  =       diagonal(k);
        }

        complex<double> exact_determinant;
        exact_determinant       =       compute_Determinant(K);
        cout << endl << "Exact log determinant is: " << setprecision(16) << exact_determinant << endl;
	//
	// cout << endl << "Exact matrix matrix product..." << endl;
	// start			=	clock();
	// MatrixXcd bExact	=	B*x;
	// end				=	clock();
	// cout << endl << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;
	//
	// cout << endl << (bExact-b).cwiseAbs().maxCoeff() << endl;
}
