#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <iomanip>
#include "HODLR_Tree.hpp"
#include "HODLR_Matrix.hpp"
#include "KDTree.hpp"

int dim = 2;

class myHODLR_Matrix : public HODLR_Matrix 
{
    private:
        int N;
        Eigen::MatrixXd Theta;

    public:

        myHODLR_Matrix(int N):HODLR_Matrix(N) 
        {
            this->N = N;
            Theta = Eigen::MatrixXd::Random(N,dim);
            get_KDTree_Sorted(Theta, 0);
        };

        double get_Matrix_Entry(int j, int k) 
        {
            // return get_Gaussian_Kernel_Entry(j,k);
            
            // double a =   2.0/pow(N,1.0/dim);
            // return get_RPY(j,k,a);
            
            return get_Biharmonic(j,k);
        }

        double get_Gaussian_Kernel_Entry(int j, int k)
        {
            if(j == k)
            {
                return 2.0;
            }

            else 
            {
                double r = 0.0;
                
                for(int i = 0; i<dim; i++) 
                {
                    r += (Theta(j, i) - Theta(k, i)) * (Theta(j, i) - Theta(k, i));
                }
                
                return exp(-r);
            }
        }

        double get_Biharmonic(int j, int k) 
        {
            if(j == k) 
            {
                return (2.0 * N);
            }
            else 
            {
                double r2 = (Theta.row(j)-Theta.row(k)).squaredNorm();
                return (0.5 * r2 * log(r2));
            }
        }

        double get_RPY(int j, int k, double a) 
        {
            if (j == k) 
            {
                return (1.0 / a);
            }
            
            else 
            {
                Eigen::VectorXd R = (Theta.row(j)-Theta.row(k));
                double RRT        = R(j%dim)*R(k%dim);
                double r          = R.norm();
                
                if (r<2.0*a) 
                {
                    double temp = 0.09375 / a * RRT / r;
                    
                    if (j%dim == k%dim) 
                    {
                        return (1.0 / a * (1 - (0.28125 * r / a) + temp));
                    }
                    
                    else 
                    {
                        return (temp / a);
                    }
                }
                
                else 
                {
                    double r2   = r*r;
                    double r4   = r2*r2;
                    double temp = RRT/r/r - 2.0*a*a/r4*RRT;
                    
                    if (j%dim==k%dim) 
                    {
                        return 0.75*(1.0+2.0/3.0*a*a/r2 + temp)/r;
                    }
                    
                    else 
                    {
                        return 0.75*temp/r;
                    }
                }
            }
        }

        ~myHODLR_Matrix() {};
};

int main(int argc, char* argv[]) 
{
    int N             = atoi(argv[1]);
    myHODLR_Matrix* A = new myHODLR_Matrix(N);
    int M             = atoi(argv[2]);
    int nLevels       = log(N/M)/log(2);
    std::cout << "\nSize of the matrix is: " << N << "\n";
    std::cout << "\nNumber of levels in the tree: " << nLevels << "\n";
    double tolerance        =   pow(10.0,-atoi(argv[3]));
    double start, end;

    double CPS  =   1.0;

    /****************/
    /*              */
    /*  Create Tree */
    /*              */
    /****************/

    std::cout << "\nCreate the tree...\n";
    start   =   omp_get_wtime();
    HODLR_Tree* myMatrix    =   new HODLR_Tree(nLevels, tolerance, A);
    end     =   omp_get_wtime();
    std::cout << "\nTree created.\n";
    double time_create_tree =   (end-start)/CPS;

    /****************************/
    /*                          */
    /*  Fast Symmetric Method   */
    /*                          */
    /****************************/

    std::cout << "\n************************************\n";
    std::cout << "\nFAST SYMMETRIC FACTORIZATION METHOD:\n";
    std::cout << "\n************************************\n";

    /********************************/
    /*                              */
    /*  Assemble Symmetric Matrix   */
    /*                              */
    /********************************/

    std::cout << "\nAssemble the matrix...\n";
    myMatrix->assembleSymmetricTree();
    end     =   omp_get_wtime();
    std::cout << "\nMatrix assembled.\n";
    double time_symmetric_assemble              =   (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Factorize Symmetric Matrix  */
    /*                              */
    /********************************/

    std::cout << "\nFactorize the matrix...\n";
    start   =   omp_get_wtime();
    myMatrix->symmetric_factorize();
    end     =   omp_get_wtime();
    std::cout << "\nFactorization done.\n";
    double time_symmetric_factorize             =   (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Determinant Symmetric Matrix*/
    /*                              */
    /********************************/

    std::cout << "\nDeterminant computation...\n";
    start   =   omp_get_wtime();
    double det  = myMatrix->symmetric_Determinant();
    end     =   omp_get_wtime();
    std::cout << "\nDeterminant computed.\n";
    double time_symmetric_determinant           =   (end-start)/CPS;

    /****************************/
    /*                          */
    /*  Solve Symmetric Matrix  */
    /*                          */
    /****************************/

    Eigen::MatrixXd b   =   Eigen::MatrixXd::Random(N,1);
    Eigen::MatrixXd xSymFast;
    std::cout << "\nSolve linear system...\n";
    start   =   omp_get_wtime();
    xSymFast=   myMatrix->symmetric_Solve(b);
    end     =   omp_get_wtime();
    std::cout << "\nSolving done.\n";
    double time_symmetric_solve                 =   (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Symmetric Matrix Product    */
    /*                              */
    /********************************/

    Eigen::MatrixXd bSymFast;
    std::cout << "\nMatrix matrix product...\n";
    start   =   omp_get_wtime();
    myMatrix->matmat_Symmetric_Product(xSymFast,bSymFast);
    end     =   omp_get_wtime();
    std::cout << "\nMatrix matrix product done.\n";
    double time_symmetric_matrix_matrix_product =   (end-start)/CPS;
    std::cout << "\nRelative error in residual when solving the linear system is: " << (bSymFast-b).norm()/b.norm() << "\n";

    /********************************/
    /*                              */
    /*  Symmetric Factor Product    */
    /*                              */
    /********************************/

    Eigen::MatrixXd bSymFacMatProd;
    std::cout << "\nSymmetric factor product...\n";
    start           =   omp_get_wtime();
    bSymFacMatProd  =   myMatrix->symmetric_Factor_Product(xSymFast);
    end             =   omp_get_wtime();
    std::cout << "\nSymmetric factor product done.\n";
    double time_symmetric_Factor_Product        =   (end-start)/CPS;

    /****************************************/
    /*                                      */
    /*  Symmetric Factor Transpose Product  */
    /*                                      */
    /****************************************/

    Eigen::MatrixXd bSymFacTranMatProd;
    std::cout << "\nSymmetric factor product...\n";
    start               =   omp_get_wtime();
    bSymFacTranMatProd  =   myMatrix->symmetric_Factor_Transpose_Product(xSymFast);
    end                 =   omp_get_wtime();
    std::cout << "\nSymmetric factor tranpose product done.\n";
    double time_symmetric_Factor_Tranpose_Product=  (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Fast Non-Symmetric Method   */
    /*                              */
    /********************************/
    std::cout << "\n****************************************\n";
    std::cout << "\nFAST NON-SYMMETRIC FACTORIZATION METHOD:\n";
    std::cout << "\n****************************************\n";

    /************************************/
    /*                                  */
    /*  Assemble Non-Symmetric Matrix   */
    /*                                  */
    /************************************/

    std::cout << "\nAssemble the matrix...\n";
    myMatrix->assemble_Tree();
    end     =   omp_get_wtime();
    std::cout << "\nMatrix assembled.\n";
    double time_assemble                =   (end-start)/CPS;

    /************************************/
    /*                                  */
    /*  Factorize Non-Symmetric Matrix  */
    /*                                  */
    /************************************/

    std::cout << "\nFactorize the matrix...\n";
    start   =   omp_get_wtime();
    myMatrix->factorize();
    end     =   omp_get_wtime();
    std::cout << "\nFactorization done.\n";
    double time_factorize               =   (end-start)/CPS;

    /************************************/
    /*                                  */
    /*  Determinant Non-Symmetric Matrix*/
    /*                                  */
    /************************************/

    std::cout << "\nDeterminant computation...\n";
    start   =   omp_get_wtime();
    double detNonSym    = myMatrix->determinant();
    end     =   omp_get_wtime();
    std::cout << "\nDeterminant computed.\n";
    double time_determinant             =   (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Solve Non-Symmetric Matrix  */
    /*                              */
    /********************************/

    Eigen::MatrixXd xFast;
    std::cout << "\nSolve linear system...\n";
    start   =   omp_get_wtime();
    xFast   =   myMatrix->solve(b);
    end     =   omp_get_wtime();
    std::cout << "\nSolving done.\n";
    double time_solve                   =   (end-start)/CPS;

    /********************************/
    /*                              */
    /*  Non-Symmetric Matrix Product*/
    /*                              */
    /********************************/

    Eigen::MatrixXd bFast;
    std::cout << "\nMatrix matrix product...\n";
    start   =   omp_get_wtime();
    myMatrix->matmat_Product(xFast,bFast);
    end     =   omp_get_wtime();
    std::cout << "\nMatrix matrix product done.\n";
    double time_matrix_matrix_product   =   (end-start)/CPS;
    std::cout << "\nRelative error in residual when solving the linear system: " << (bFast-b).norm()/b.norm() << "\n";

    std::cout << "\nRelative error in residual comparing the symmetric and non-symmetric solution: " << (bFast-bSymFast).norm()/bFast.norm() << "\n";
    std::cout << "\nRelative error in the log-determinant using the symmetric and non-symmetric method: " << fabs(detNonSym-det)/fabs(detNonSym) << "\n";

    /************/
    /*          */
    /*  Summary */
    /*          */
    /************/
    std::cout << "\n*********************\n";
    std::cout << "\nSUMMARY OF TIME TAKEN\n";
    std::cout << "\n*********************\n";
    std::cout << std::setw(20);
    std::cout << "\n\tCreate tree: " << time_create_tree << "\n";
    std::cout << "\nFOR SYMMETRIC FACTORIZATION: \n";
    std::cout << "\n\tAssemble matrix: " << time_symmetric_assemble << "\n";
    std::cout << "\n\tFactorize matrix: " << time_symmetric_factorize << "\n";
    std::cout << "\n\tComputing determinant: " << time_symmetric_determinant << "\n";
    std::cout << "\n\tSolving linear system: " << time_symmetric_solve << "\n";
    std::cout << "\n\tSymmetric matrix vector product: " << time_symmetric_matrix_matrix_product << "\n";
    std::cout << "\n\tSymmetric factor product: " << time_symmetric_Factor_Product << "\n";
    std::cout << "\n\tSymmetric factor transpose product: " << time_symmetric_Factor_Tranpose_Product << "\n";
    std::cout << "\nFOR NON-SYMMETRIC FACTORIZATION: \n";
    std::cout << "\n\tAssemble matrix: " << time_assemble << "\n";
    std::cout << "\n\tFactorize matrix: " << time_factorize << "\n";
    std::cout << "\n\tComputing determinant: " << time_determinant << "\n";
    std::cout << "\n\tSolving linear system: " << time_solve << "\n";
    std::cout << "\n\tMatrix vector product: " << time_matrix_matrix_product << "\n";
}
