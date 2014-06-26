//
//  HODLR_diagnostics.cpp
//  
//
//  Created by Sivaram Ambikasaran on 2/7/14.
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
        
	unsigned N	=	10000;
	unsigned nRhs	=	1;
	unsigned nLeaf	=	50;
	double tolerance=	1e-16;
        
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
	VectorXd diagonal	=	2.0*VectorXd::Ones(N);
	A->assemble_Matrix(diagonal, tolerance);
	end		=	clock();
	cout << "Time taken is: " << double(end-start)/double(CLOCKS_PER_SEC)<< endl;

        A->diagnostics();

        //      Displays the off-diagonal rank at all levels and all nodes.
        A->display_all_Ranks();

        //      Displays the off-diagonal rank of the '4'th matrix at the '3'rd level. For instance, A->display_Rank(0,0) will display the rank of the largest two off-diagonal blocks.
        A->display_Rank(3,4);

        //      Displays an estimate of total number of calls to get_Matrix_Entry.
        A->total_Calls_To_Get_Matrix_Entry();

        //      Displays an estimate of number of calls to get_Matrix_Entry of the '9'th matrix at the '4'th level.
        A->total_Calls_To_Get_Matrix_Entry(4,9);
}