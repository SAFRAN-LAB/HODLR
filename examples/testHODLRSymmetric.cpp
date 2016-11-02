#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

int dim  = 1;
double tol = 1e-13;

class myHODLR_Matrix : public HODLR_Matrix {
	private:
		//Eigen::VectorXd x;
		Eigen::MatrixXd Theta;
	public:
		myHODLR_Matrix(int N) : HODLR_Matrix(dim*N) {
			//x	=	Eigen::VectorXd::Random(N);
			//std::sort(x.data(),x.data()+x.size());
			Theta = Eigen::MatrixXd::Random(N,dim);
			get_KDTree_Sorted(Theta, 0);
		};

		double get_Matrix_Entry(int j, int k) {
 //           return get_Biharmonic_Kernel_Entry(j,k);
            return get_Gaussian_Kernel_Entry(j,k);
//		     return get_RPY_Tensor_Entry(j,k);
		}

        double get_Biharmonic_Kernel_Entry(int j, int k){
            if (j == k){
                return 1;
            } else {
                 Eigen::VectorXd X = Theta.row(j/dim) - Theta.row(k/dim);
                 double r = 0.0;
                 for(int i = 0; i<dim; i++){
                    r += X(i)*X(i);
                 }
                 double ex = 0.5 * r * log(r);
               if(j%dim == k%dim){
                    return 1 + ex;
                } else {
                    return ex;
                }
            }
        }

		double get_Gaussian_Kernel_Entry(int j, int k){
		     if(j == k){
		        return 20;
		     } else{
	             Eigen::VectorXd X = Theta.row(j/dim) - Theta.row(k/dim);
                 double r = 0.0;
                 for(int i = 0; i<dim; i++){
                    r += X(i)*X(i);
                 }
                 double ex = exp(-r);
                 if (j%dim == k%dim){
                    return (1 + ex);
                 } else {
                    return ex;
                 }
		     }
		}

		double get_RPY_Tensor_Entry(int j, int k){
             Eigen::VectorXd X = Theta.row(j/dim) - Theta.row(k/dim);
             double R = X(j%dim) * X(k%dim);
             double r = 0.0;
             for(int i = 0; i<dim; i++){
                r += X(i)*X(i);
             }
             r = sqrt(r);
             double a = pow(1/1500, 1/dim);
                if(r < 2*a){
                    if(j == k){
                        return 4/(3*a);
                    } else {
                        double A = 3/32 * R/a * R/r;
                        if(j%dim == k%dim){
                           return 4/3 * 1/a * (1 - 9/32 * r/a + A);
                        } else{
                            return 4/3 * A/a;
                        }
                    }
                } else {
                    if(j == k){
                        return 1/r;
                    } else {
                        double B = (R/r * R/r) - (2/3 * a/r *a/r) * (3 * R/r * R/r);
                        if(j%dim == k%dim){
                            return  1/r * (1 + 2/3 * a/r * a/r + B);
                        } else {
                            return B/r;
                        }
                    }
                }
		}
		~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) {

	int N					=	atoi(argv[1]);
	myHODLR_Matrix* A		=	new myHODLR_Matrix(N);
	int nLevels				=	log((dim*N)/200)/log(2);
	std::cout<<nLevels<<"\n";
	double tolerance		=	tol;
	double start, end;
	double CPS	=	1.0;

	std::cout << "\nFast method...\n";
	start	=	omp_get_wtime();
	HODLR_Tree* myMatrix	=	new HODLR_Tree(nLevels, tolerance, A);
	Eigen::MatrixXd x		=	Eigen::MatrixXd::Random(dim*N,1);
	Eigen::MatrixXd bFast;
	myMatrix->assembleSymmetricTree();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken for assembling the matrix in HODLR form is: " << (end-start)/CPS << "\n";

 //   std::cout<<" rank: "<<myMatrix->tree[0][0]->sym_rank<<"\n";

	Eigen::MatrixXd xFast;
	start	=	omp_get_wtime();
	myMatrix->symmetric_factorize();
	end		=	omp_get_wtime();
	std::cout << "\nTime taken to factorize is: " << (end-start)/CPS << std::endl;




    Eigen::MatrixXd M = A->get_Matrix(0,0,N,N);
    Eigen::LLT<Eigen::MatrixXd> P;
    P.compute(M);
    double det = 0.0;
    for(int i=0; i<P.matrixL().rows(); ++i){
        det += log(P.matrixL()(i,i));
    }
       start	=	omp_get_wtime();
       double det1 = myMatrix->symmetric_Determinant();
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_determinant: " << (end-start)/CPS << "\n";
       std::cout<<"fast determinant: "<<det1<<" Usual det: "<<2*det<<"\n";

       /* start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Factor_Solve(x);
        end		=	omp_get_wtime();
        std::cout << "\nTime taken for symmetric_Inverse: " << (end-start)/CPS << "\n";
        std::cout<<"bfast norm "<<(x - Rc*bFast).norm()<<"\n";

        start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Factor_Transpose_Solve(x);
        end		=	omp_get_wtime();
        std::cout << "\nTime taken for symmetric_Transpose_Inverse: " << (end-start)/CPS << "\n";
       std::cout<<"bfast norm "<<(x - Rc.transpose()*bFast).norm()<<"\n";*/

        /*start	=	omp_get_wtime();
        bFast = myMatrix->symmetric_Solve(x);
        end		=	omp_get_wtime();
        std::cout << "\nTime taken for symmetric_Full_Inverse: " << (end-start)/CPS << "\n";
        std::cout<<"bfast norm "<<(x - M*bFast).norm()<<"\n";*/

       /*start	=	omp_get_wtime();
       bFast = myMatrix->symmetric_Factor_Product(x);
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_Product: " << (end-start)/CPS << "\n";
       std::cout<<"bfast norm "<<(Rc*x - bFast).norm()<<"\n";

       start	=	omp_get_wtime();
       bFast = myMatrix->symmetric_Factor_Transpose_Product(x);
       end		=	omp_get_wtime();
       std::cout << "\nTime taken for symmetric_Transpose_Product: " << (end-start)/CPS << "\n";
       std::cout<<"bfast norm "<<(Rc.transpose()*x - bFast).norm()<<"\n";*/

       start = omp_get_wtime();
       bFast = myMatrix->build_Symmetric_Matrix_Factor();
       end = omp_get_wtime();
       std::cout<<"\nTime taken for building matrix: "<<(end-start)/CPS<<"\n";
       std::cout<<"bfast norm "<<(M - bFast*bFast.transpose()).norm()<<"\n";


    }