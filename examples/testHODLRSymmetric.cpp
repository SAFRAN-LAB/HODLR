#include <iostream>
#include <Eigen/Dense>
#include <cmath>
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
				return 100.0;
			}
			else {
//				return 1.0/(1.0+((x(j)-x(k))*(x(j)-x(k))));
				return exp(-(x(j)-x(k))*(x(j)-x(k)));
//                return exp(-fabs(x(j)-x(k)));
			}
		}
		~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) {
	int N					=	atoi(argv[1]);
	myHODLR_Matrix* A		=	new myHODLR_Matrix(N);
//	int nLevels				=	log(N/200)/log(2);
	int nLevels				=	log(N/256)/log(2);
	double tolerance		=	1e-15;



	double start, end;
	// double CPS	=	CLOCKS_PER_SEC;
	double CPS	=	1.0;

	std::cout << "\nFast method...\n";
	// start	=	clock();
	start	=	omp_get_wtime();
	HODLR_Tree* myMatrix	=	new HODLR_Tree(nLevels, tolerance, A);
	Eigen::MatrixXd x		=	Eigen::MatrixXd::Random(N,1);
	Eigen::MatrixXd bFast;
	myMatrix->assembleSymmetricTree();
	// end		=	clock();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for assembling the matrix in HODLR form is: " << (end-start)/CPS << "\n";

	Eigen::MatrixXd xFast;
	// start	=	clock();
	start	=	omp_get_wtime();
	myMatrix->symmetric_factorize();
	// end	=	clock();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken to factorize is: " << (end-start)/CPS << std::endl;


    Eigen::MatrixXd M;
  Eigen::MatrixXd Wc;
    Eigen::MatrixXd Rc = Eigen::MatrixXd::Identity(N,N);
    for(int l=0;l<nLevels;++l){
        Wc = Eigen::MatrixXd::Zero(N,N);
        for(int k=0; k<myMatrix->nodesInLevel[l]; ++k){
                   int r = myMatrix->tree[l][k]->sym_rank;
                   int n0 = myMatrix->tree[l][k]->Q[0].rows();
                   int n1 = myMatrix->tree[l][k]->Q[1].rows();
                   Eigen::MatrixXd T = Eigen::MatrixXd::Identity(n0+n1,n0+n1);
                   Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(n0+n1,2*r);
                   temp.block(0,0,n0,r) = myMatrix->tree[l][k]->Qfactor[0];
                   temp.block(n0,r,n1,r) = myMatrix->tree[l][k]->Qfactor[1];
                   T += temp*myMatrix->tree[l][k]->X*temp.transpose();
                   Wc.block(myMatrix->tree[l][k]->nStart, myMatrix->tree[l][k]->nStart, myMatrix->tree[l][k]->nSize, myMatrix->tree[l][k]->nSize) = T;

        }
 //     Rc = Wc*Rc*Wc.transpose();
        Rc = Wc*Rc;
    }

    Wc = Eigen::MatrixXd::Zero(N,N);


    int l=nLevels;
    	for(int k=0; k<myMatrix->nodesInLevel[l]; ++k){
        Wc.block(myMatrix->tree[l-1][k/2]->cStart[k%2],myMatrix->tree[l-1][k/2]->cStart[k%2],myMatrix->tree[l-1][k/2]->cSize[k%2],myMatrix->tree[l-1][k/2]->cSize[k%2]) = myMatrix->tree[l][k]->llt.matrixL();
    }


//    Rc = Wc*Rc*Wc.transpose();
    Rc = Wc*Rc;

    Eigen::MatrixXd Rc1 = Rc*Rc.transpose();
	 M = A->get_Matrix(0,0,N,N);
//	Eigen::LLT<Eigen::MatrixXd> mllt;
//    mllt.compute(M);
//    Eigen::MatrixXd L = mllt.matrixL();

//    Eigen::MatrixXd P = Rc-M;


	/*double nrm = ((double)(Rc1 - M).norm());
    	std::cout<<"Error norm: "<<std::endl<<nrm<<std::endl;*/

   // double det = 0.0;

    /*for(int i= 0; i<L.rows(); ++i){
        det += log(L(i,i))/log(2);
    }*/
    	/*std::cout<<"Determinant "<<myMatrix->symmetric_determinant()<<" Det "<< 2*det<<"\n";*/

       start	=	omp_get_wtime();
       /*double det =*/ myMatrix->symmetric_Determinant();
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_determinant: " << (end-start)/CPS << "\n";

        start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Inverse(x);
        end		=	omp_get_wtime();
        std::cout<<"bfast norm "<<(x - Rc*bFast).norm()<<"\n";
        std::cout << "\nTime taken for symmetric_Inverse: " << (end-start)/CPS << "\n";

        start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Transpose_Inverse(x);
        end		=	omp_get_wtime();
        std::cout << "\nTime taken for symmetric_Transpose_Inverse: " << (end-start)/CPS << "\n";
        std::cout<<"bfast norm "<<(x - Rc.transpose()*bFast).norm()<<"\n";

        start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Full_Inverse(x);
        end		=	omp_get_wtime();
        std::cout << "\nTime taken for symmetric_Full_Inverse: " << (end-start)/CPS << "\n";
        std::cout<<"bfast norm "<<(x - M*bFast).norm()<<"\n";

       start	=	omp_get_wtime();
       bFast = myMatrix->symmetric_Product(x);
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_Product: " << (end-start)/CPS << "\n";
       std::cout<<"bfast norm "<<(Rc*x - bFast).norm()<<"\n";

       start	=	omp_get_wtime();
       bFast = myMatrix->symmetric_Transpose_Product(x);
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_Transpose_Product: " << (end-start)/CPS << "\n";
       std::cout<<"bfast norm "<<(Rc.transpose()*x - bFast).norm()<<"\n";

}
