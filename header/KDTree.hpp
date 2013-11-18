//
//  KDTree.hpp
//  
//
//  Created by Sivaram Ambikasaran on 10/22/13.
//
//

#ifndef __KDTREE_HPP__
#define __KDTREE_HPP__

#include"iostream"
#include"Eigen/Dense"

using namespace std;
using namespace Eigen;

void mergeSortedLists(MatrixXd& list1, MatrixXd& list2, unsigned index, MatrixXd& finalList);

void mergeSort(MatrixXd& locations, unsigned index);

void get_KDTree_Sorted(MatrixXd& locations, unsigned index);

#endif /* defined(__KDTREE_HPP__) */
