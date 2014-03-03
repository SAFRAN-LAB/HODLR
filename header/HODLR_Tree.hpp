//
//  HODLR_Tree.hpp
//  
//
//  Created by Sivaram Ambikasaran on 3/2/14.
//
//

#ifndef __HODLR_Tree_hpp__
#define __HODLR_Tree_hpp__

#include "helper_Functions.hpp"
#include <vector>
#include <cmath>
#include <Eigen/Dense>

using namespace Eigen;

template <typename MatrixType>
class HODLR_Tree {
public:
        /**     */
        friend class HODLR_Node<MatrixType>;
        MatrixType matrix;
        vector< vector<HODLR_Node<MatrixType> > > tree;
        int N;
        int nLeaf;
        double lowRankTolerance;
        double determinant;
        int nLevels;
        int nNodes;

        /*!
         Constructor for HODLR tree.
         */
        HODLR_Tree(MatrixType* matrix, int N, int nLeaf) {
                this->matrix    =       matrix;
                this->N         =       N;
                this->nLeaf     =       nLeaf;
                this->nLevels   =       log2(N/nLeaf);
                this->nNodes    =       power(2,nLevels+1)-1;
                create_Tree();
        }
        
        /*!
         Creates the tree.
         */
        void create_Tree() {
                
        }
};

#endif /*(__HODLR_tree_hpp__)*/