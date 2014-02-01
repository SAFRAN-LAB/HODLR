/*!
 \brief This function is for the construction of KDTree.
 
 \note
 
 \author $Sivaram Ambikasaran$
 
 \version
 
 \date $Date: November 14th, 2013$
 */

#include "KDTree.hpp"
#include<Eigen/Dense>

using namespace Eigen;

void mergeSortedLists(MatrixXd& list1, MatrixXd& list2, unsigned index, MatrixXd& finalList) {
	unsigned N1	=	list1.rows();
	unsigned N2	=	list2.rows();
	unsigned j1	=	0;
	unsigned j2	=	0;
	unsigned j	=	0;
	while (j1<N1 && j2 <N2) {
		if (list1(j1,index) < list2(j2,index)) {
			finalList.row(j)	=	list1.row(j1);
			++j1;
		}
		else {
			finalList.row(j)	=	list2.row(j2);
			++j2;
		}
		++j;
	}
	while (j1<N1) {
		finalList.row(j)	=	list1.row(j1);
		++j1;
		++j;
	}
	while (j2<N2) {
		finalList.row(j)	=	list2.row(j2);
		++j2;
		++j;
	}
}

void mergeSort(MatrixXd& locations, unsigned index) {
	unsigned N		=	locations.rows();
	if (N==1) {
		return;
	}
	else {
		///	Number of points in the left cluster.
		unsigned Nleft		=	N/2;
		
		///	Number of points in the right cluster.
		unsigned Nright		=	N-Nleft;
		
		///	Dimension of the space.
		unsigned nDimensions	=	locations.cols();
		
		///	Left locations.
		MatrixXd leftLocations	=	locations.block(0,0,Nleft,nDimensions);
		
		///	Right locations.
		MatrixXd rightLocations	=	locations.block(Nleft,0,Nright,nDimensions);
		
		///	Mergesort for the left.
		mergeSort(leftLocations, index);
		
		///	Mergesort for the right.
		mergeSort(rightLocations, index);
		
		///	Merge the sorted left and right lists.
		mergeSortedLists(leftLocations, rightLocations, index, locations);
	}
}

void get_KDTree_Sorted(MatrixXd& locations, unsigned index) {
	///	Get the total number of points.
	unsigned N		=	locations.rows();
	if (N==1) {
		return;
	}
	else {
		///	Number of points in the left cluster.
		unsigned Nleft		=	N/2;

		///	Number of points in the right cluster.
		unsigned Nright		=	N-Nleft;

		///	Dimension of the space.
		unsigned nDimensions	=	locations.cols();

		///	Merge sort on the input locations based on the coordinate index%2.
		mergeSort(locations, index%nDimensions);

		///	Obtain the left and right locations.
		MatrixXd leftLocations	=	locations.block(0,0,Nleft,nDimensions);
		MatrixXd rightLocations	=	locations.block(Nleft,0,Nright,nDimensions);

		///	Sort the left and right locations based on a KDTree.
		get_KDTree_Sorted(leftLocations, index+1);
		get_KDTree_Sorted(rightLocations, index+1);

		///	Output the locations.
		locations.block(0,0,Nleft,nDimensions)		=	leftLocations;
		locations.block(Nleft,0,Nright,nDimensions)	=	rightLocations;
	}
}