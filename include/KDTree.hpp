#ifndef __KDTree__
#define __KDTree__

#include "Eigen/Dense"
#include "HODLR_Matrix.hpp"

void mergeSortedLists(Mat& list1, Mat& list2, unsigned index, Mat& final_list);
void mergeSort(Mat& locations, unsigned index);
void getKDTreeSorted(Mat& locations, unsigned index);

#endif /*__KDTree__*/
