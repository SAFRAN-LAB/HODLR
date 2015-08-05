#include<iostream>
#include<Eigen/Dense>

using namespace std;

int main() {
	int m = 100;
	Eigen::MatrixXd A	=	Eigen::MatrixXd::Zero(m,m);
	std::cout << A << std::endl;
}