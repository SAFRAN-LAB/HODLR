#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

int dim     = 2;
double tol  = 1e-12;

class myHODLR_Matrix : public HODLR_Matrix {
	private:
		int N;
		Eigen::MatrixXd Theta;
	public:
		myHODLR_Matrix(int N) : HODLR_Matrix(N) {
			this->N = N;
			Theta = Eigen::MatrixXd::Random(N,dim);
			get_KDTree_Sorted(Theta, 0);
		};

		double get_Matrix_Entry(int j, int k) {
			return get_Gaussian_Kernel_Entry(j,k);
		}

		double get_Gaussian_Kernel_Entry(int j, int k){
			if(j == k){
				return 2.0;
			}
			else {
				double r = 0.0;
				for(int i = 0; i<dim; i++) {
					r += (Theta(j, i) - Theta(k, i)) * (Theta(j, i) - Theta(k, i));
				}
				return exp(-r);
			}
		}

		~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) {

	int N					=	atoi(argv[1]);
	myHODLR_Matrix* A		=	new myHODLR_Matrix(N);
	int nLevels				=	log(N/200)/log(2);
	std::cout<<nLevels<<"\n";
	double tolerance		=	tol;
	double start, end;
	double CPS	=	1.0;

	std::cout << "\nFast method...\n";
	start	=	omp_get_wtime();
	HODLR_Tree* myMatrix	=	new HODLR_Tree(nLevels, tolerance, A);
	Eigen::MatrixXd x		=	Eigen::MatrixXd::Random(N,1);
	Eigen::MatrixXd bFast, bFast1;
	myMatrix->assembleSymmetricTree();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for assembling the matrix in HODLR form is: " << (end-start)/CPS << "\n";

	Eigen::MatrixXd xFast, xFast1;
	start	=	omp_get_wtime();
	myMatrix->symmetric_factorize();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken to factorize is: " << (end-start)/CPS << std::endl;



	// Don't use it for large N. Mainly used for testing purposes.
	Eigen::MatrixXd M = A->get_Matrix(0,0,N,N);
	Eigen::LLT<Eigen::MatrixXd> P;
	P.compute(M);
	double det = 0.0;
	for(int i=0; i<P.matrixL().rows(); ++i){
		det += log(P.matrixL()(i,i));
	}
	det = 2*det;



	start	=	omp_get_wtime();
	double det1 = myMatrix->symmetric_Determinant();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for symmetric_determinant: " << (end-start)/CPS << "\n";
	std::cout<<"fast determinant: "<<det1<<" slow: "<<det<<"\n";

	start = omp_get_wtime();
	bFast = myMatrix->symmetric_Solve(x);
	end		=	omp_get_wtime();
	xFast = M*bFast;
	std::cout<<"Error norm: "<<(xFast-x).norm()/x.norm()<<"\n";

	start	=	omp_get_wtime();
	bFast = myMatrix->symmetric_Factor_Product(x);
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for symmetric_Factor_Product: " << (end-start)/CPS << "\n";
	//   std::cout<<"bfast norm "<<(Rc*x - bFast).norm()<<"\n";

	start	=	omp_get_wtime();
	bFast = myMatrix->symmetric_Factor_Transpose_Product(x);
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for symmetric_Factor_Transpose_Product: " << (end-start)/CPS << "\n";
	//   std::cout<<"bfast norm "<<(Rc.transpose()*x - bFast).norm()<<"\n";

	/*start = omp_get_wtime();
	bFast = myMatrix->build_Symmetric_Matrix_Factor();
	end = omp_get_wtime();
	std::cout<<"\nTime taken for building matrix: "<<(end-start)/CPS<<"\n";
	std::cout<<"bfast norm "<<(M - bFast*bFast.transpose()).norm()<<"\n";*/

	start = omp_get_wtime();
	myMatrix->assembleTree();
	end		=	omp_get_wtime();
	std::cout << "\nHODLR non-symmetric method Time taken to assemble is: " << (end-start)/CPS << std::endl;

	start	=	omp_get_wtime();
	myMatrix->factorize();
	end		=	omp_get_wtime();
	std::cout << "\nHODLR non-symmetric Time taken to factorize is: " << (end-start)/CPS << std::endl;

	start = omp_get_wtime();
	bFast1 = myMatrix->solve(x);
	end		=	omp_get_wtime();
	xFast1 = M*bFast1;
	std::cout<<"HODLR non-symmetric method Error norm: "<<(xFast1-x).norm()/x.norm()<<"\n";
}
