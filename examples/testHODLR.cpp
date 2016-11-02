#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"

class myHODLR_Matrix : public HODLR_Matrix {
private:
	Eigen::VectorXd x;
public:
	myHODLR_Matrix(int N) : HODLR_Matrix(N) {
		x	=	Eigen::VectorXd::Random(N);
		std::sort(x.data(),x.data()+x.size());
	};
	double get_Matrix_Entry(int j, int k) {
		if(j==k) {
			return 10;
		}
		else {
			//return 1.0/(1.0+((x(j)-x(k))*(x(j)-x(k))));
//			return sqrt(1.0 + ((x(j)-x(k))*(x(j)-x(k))));
			// return exp(-fabs(x(j)-x(k)));
			return exp(-(x(j)-x(k))*(x(j)-x(k)));
		}
	}
	~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) {
	int N					=	atoi(argv[1]);
	myHODLR_Matrix* A		=	new myHODLR_Matrix(N);
	int nLevels				=	log(N/200)/log(2);
	double tolerance		=	1e-13;

	double start, end;
	// double CPS	=	CLOCKS_PER_SEC;
	double CPS	=	1.0;

	std::cout << "\nFast method...\n";
	// start	=	clock();
	start	=	omp_get_wtime();
	HODLR_Tree* myMatrix	=	new HODLR_Tree(nLevels, tolerance, A);
	Eigen::MatrixXd x		=	Eigen::MatrixXd::Random(N,1);
	Eigen::MatrixXd bFast;
	myMatrix->assembleTree();
	// end		=	clock();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for assembling the matrix in HODLR form is: " << (end-start)/CPS << "\n";

	// start	=	clock();
	start	=	omp_get_wtime();
	myMatrix->matmat_Product(x, bFast);
	// end		=	clock();
	end		=	omp_get_wtime();
	std::cout << "\nTime for fast matrix-vector product is: " << (end-start)/CPS << "\n";


	// std::cout << "\nExact method...\n";
	// // start	=	clock();
	// start	=	omp_get_wtime();
	// Eigen::MatrixXd B		=	A->get_Matrix(0,0,N,N);
	// // end		=	clock();
	// end		=	omp_get_wtime();
	// std::cout << "\nTime taken for assembling the matrix is: " << (end-start)/CPS << "\n";

	// // start	=	clock();
	// start	=	omp_get_wtime();
	// Eigen::MatrixXd bExact	=	B*x;
	// // end	=	clock();
	// end	=	omp_get_wtime();
	//
	// std::cout << "\nTime for exact matrix-vector product is: " << (end-start)/CPS << "\n";
	//
	// std::cout << "\nError in the solution is: " << (bFast-bExact).norm()/(1.0+bExact.norm()) << "\n";

	Eigen::MatrixXd xFast;
	// start	=	clock();
	start	=	omp_get_wtime();
	myMatrix->factorize();
	// end	=	clock();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken to factorize is: " << (end-start)/CPS << "\n";

	// start	=	clock();
	start	=	omp_get_wtime();
	xFast	=	myMatrix->solve(bFast);
	// end	=	clock();
	end	=	omp_get_wtime();
	std::cout << "\nTime taken to solve is: " << (end-start)/CPS << "\n";

	// start	=	clock();
	// Eigen::MatrixXd xSolve	=	B.fullPivLu().solve(bExact);
	// end		=	clock();
	// std::cout << "\nTime taken to solve is: " << (end-start)/CPS << "\n";

	std::cout << "\nError in the solution is: " << (xFast-x).norm()/*/(1.0+x.norm())*/ << "\n";

  Eigen::MatrixXd M = A->get_Matrix(0,0,N,N);
Eigen::LLT<Eigen::MatrixXd> P;
    P.compute(M);
    double det = 0.0;
    for(int i=0; i<P.matrixL().rows(); ++i){
        det += log(P.matrixL()(i,i));
    }

	start = omp_get_wtime();
           double det1 = myMatrix->hodlr_Determinant();
           end = omp_get_wtime();
           std::cout<<"\nTime taken for computing: "<<(end-start)/CPS<<"\n";

           std::cout<<"fast determinant: "<<det1<<" Usual det: "<<2*det<<"\n";
}
