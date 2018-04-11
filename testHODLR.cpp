//version	:	1.1.1
//date		:	11-4-2018(dd-mm-yyyy)
//author	:	Vaishnavi G.

#include <cstdlib>
#include "HODLR.hpp"
#include <ctime>


class userkernel: public kernel {
public:
	#ifdef ONEOVERR
	userkernel() {
		alpha		=	-1.0;
	};
	double getInteraction(const double r1, const double r2) {
		double R	=	fabs(r1-r2);
		return 1.0/(R+0.5);
	};
	#elif LOGR
	userkernel() {
		alpha		=	1.0;
	};
	double getInteraction(const double r1, const double r2) {
		double R	=	fabs(r1-r2);
		return log(R+0.5);
	};
	#elif INVEXPR
	userkernel() {
		alpha		=	1.0;
	};
	double getInteraction(const double r1, const double r2) {
		double R	=	fabs(r1-r2);
		double alpha	=	2.0;
		return exp(-alpha*R);
	};
	#endif
	Eigen::VectorXd charges;
	~userkernel() {};
};


int main(int argc, char* argv[]) {
	double L	=	std::strtod (argv[1],NULL);	//semiLength of rod
	int nChebNodes	=	atoi(argv[2]);	//no. of chebnodes
	int nLevels	=	atoi(argv[3]);  //no. of levels
	double start, end;
	int N = pow(2,nLevels+1)*nChebNodes; //no. of charges
	Eigen::VectorXd charges;
	
	Eigen::VectorXd b = Eigen::VectorXd::Random(N); //right hand side vector

	start	=	omp_get_wtime();
	userkernel* mykernel		=	new userkernel();
	HODLR1DTree<userkernel>* Q	=	new HODLR1DTree<userkernel>(mykernel, nLevels, nChebNodes, N, b, charges, L);

	Q->set_Standard_Cheb_Nodes();

	Q->set_Standard_Charges();

	Q->createTree();

	Q->assign_Center_Location();

	Q->assign_Leaf_Charges();

	Q->assign_ChebNodes();

	
	end		=	omp_get_wtime();
	double timeCreateTree	=	(end-start);

	std::cout << std::endl << "Number of particles is: " << Q->N << std::endl;
	std::cout << std::endl << "Time taken to create the tree is: " << timeCreateTree << std::endl;


	start	=	omp_get_wtime();

	Q->assembleMatrix();

	end		=	omp_get_wtime();
	double timeAssembleMatrix=	(end-start);
	std::cout << std::endl << "Time taken to assemble the matrix is: " << timeAssembleMatrix << std::endl;

	start	=	omp_get_wtime();
	
	Q->solve();

	end		=	omp_get_wtime();
	double timeSolve			=	(end-start);
	std::cout << std::endl << "Time taken to solve is: " << timeSolve << std::endl;


	double totalTime	=	timeCreateTree+timeAssembleMatrix+timeSolve;
	
	double applyTime	=	timeAssembleMatrix+timeSolve;

	std::cout << std::endl << "Total time taken is: " << totalTime << std::endl;

	std::cout << std::endl << "Apply time taken is: " << applyTime << std::endl;

	std::cout << std::endl << "Total Speed in particles per second is: " << Q->N/totalTime << std::endl;

	std::cout << std::endl << "Apply Speed in particles per second is: " << Q->N/applyTime << std::endl;

	std::cout << std::endl << "Performing Error check..." << std::endl;

	srand(time(NULL));
	int nLine	=	rand()%Q->nLinesPerLevel[nLevels];
	std::cout << std::endl << "For Line number: " << nLine << std::endl;
	
	start	=	omp_get_wtime();
	std::cout << std::endl << "Error is: " << Q->perform_Error_Check(nLine) << std::endl;
	end		=	omp_get_wtime();
	double ErrorTime	=	(end-start);

	std::cout << std::endl << "Time taken to compute error is: " << ErrorTime << std::endl;

	std::cout << std::endl;
}
