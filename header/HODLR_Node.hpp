//
//  HODLR_Node.hpp
//  
//
//  Created by Sivaram Ambikasaran on 3/2/14.
//
//

#ifndef __HODLR_Node_hpp__
#define __HODLR_Node_hpp__

#include<Eigen/Dense>

template <typename MatrixType>
class HODLR_Node {
private:
        MatrixType matrix;
        
public:
        /**     Variables of interest for both leaf and non-leaf;       */
        int parent;                     ///<    levelBasedNodeNumber of its parent;
        int child[2];                   ///<    levelBasedNodeNumber of its children;

        int levelNumber;                ///<    Level number of the node;
        int nodeNumber;                 ///<    Node number is either '0' or '1' depending on left or right child;
        int levelBasedNodeNumber;       ///<    Level based node number;
        int nStart;                     ///<    Starting row/column index of the sub-matrix corresponding to this node;
        int nSize;                      ///<    Size (i.e., the number of rows/columns) of the sub-matrix corresponding to this node;

        MatrixXd Kself;                 ///<    Stores the self interaction at the leaf level; At the non-leaf stores the matrix [I, Vinverse[1]*Uinverse[1]; Vinverse[0]Uinverse[0], I];
        FullPivLU<MatrixXd> Kinverse;   ///<    Stores factorization of K;
        double determinant;             ///<    Stores the determinant of K, i.e., K.determinant();

        /**     Variables of interest at non-leaf;       */
        int nRank[2];                   ///<    Rank of the off-diagonal blocks for the node.
        MatrixXd U[2];                  ///<    Column basis of low-rank interaction of children at non-leaf node;
        MatrixXd V[2];                  ///<    Row basis of low-rank interaction of children at non-leaf node;

        MatrixXd Uinverse[2];           ///<    Column basis of low-rank interaction of children of inverse at non-leaf node;
        MatrixXd Vinverse[2];           ///<    Row basis of low-rank interaction of children of inverse at non-leaf node;

        /**     Variables of interest at leaf;       */
        bool isLeaf;

        /*!
         Constructor for the class.
         */
        HODLR_Node(MatrixType matrix, int levelNumber, int nodeNumber, int levelBasedNodeNumber, int nStart, int nSize) {
                this->matrix                    =       matrix;
                this->levelNumber               =       levelNumber;
                this->nodeNumber                =       nodeNumber;
                this->levelBasedNodeNumber      =       levelBasedNodeNumber;
                this->nStart                    =       nStart;
                this->nSize                     =       nSize;
                this->parent                    =       -1;
                this->child[0]                  =       -1;
                this->child[1]                  =       -1;
        }
};

#endif /*(__HODLR_Node_hpp__)*/