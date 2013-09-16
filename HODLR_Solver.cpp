/*!
 *  \copyright This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *  \author Sivaram Ambikasaran
 *  \version 3.1
 */
/*! \file HODLR_Solver.cpp
 */

#include "HODLR_Solver.hpp"

MatrixXd K;

void assemble_Matrix(const unsigned M, const unsigned N){
    /****************************/
    /*      Random Matrix       */
    /****************************/
    #ifdef RANDOM
    cout << "Matrix is a random matrix." << endl;
    K   =   MatrixXd::Random(M,N);

    /****************************/
    /*      K   =   1.0/R       */
    /****************************/
    #elif ONEOVERR
    cout << "Matrix is 1.0/R." << endl;
    VectorXd theta  =   VectorXd::Random(M);
    sort(theta.data(), theta.data()+theta.size());
    K   =   MatrixXd(M,N);
    double R;
    double scale    =   1;
    
    for (unsigned j=0; j<M; ++j) {
        for (unsigned k=0; k<N; ++k) {
            R   =   fabs(theta(j)-theta(k));
            if (R==0) {
                K(j,k)   =   scale;
            }
            else {
                K(j,k)   =   1.0/R;
            }
        }
    }
    
    /****************************/
    /*      K   =   exp(-R)     */
    /****************************/
    #elif EXPONENTIAL
    cout << "Matrix is exp(-R)." << endl;
    VectorXd theta  =   VectorXd::Random(M);
    sort(theta.data(), theta.data()+theta.size());
    K   =   MatrixXd(M,N);
    double R;
    double scale    =   1;
    for (unsigned j=0; j<M; ++j) {
        for (unsigned k=0; k<N; ++k) {
            R   =   fabs(theta(j)-theta(k));
            if (R==0) {
                K(j,k)   =   scale;
            }
            else {
                K(j,k)   =   log(R);
            }
        }
    }
    
    /****************************/
    /*      K   =   log(R)      */
    /****************************/
    #elif LOGARITHM
    cout << "Matrix is log(R)." << endl;
    VectorXd theta  =   VectorXd::Random(M);
    sort(theta.data(), theta.data()+theta.size());
    K   =   MatrixXd(M,N);
    double R;
    double scale    =   1;
    for (unsigned j=0; j<M; ++j) {
        for (unsigned k=0; k<N; ++k) {
            R       =   fabs(theta(j)-theta(k));
            K(j,k)  =   exp(-scale*R);
        }
    }
    
    #elif GAUSSIAN
    /****************************/
    /*      K   =   exp(-R^2)   */
    /****************************/
    cout << "Matrix is exp(-R^2)." << endl;
    VectorXd theta  =   VectorXd::Random(M);
    sort(theta.data(), theta.data()+theta.size());
    K   =   MatrixXd(M,N);
    double R;
    double scale    =   1;
    for (unsigned j=0; j<M; ++j) {
        for (unsigned k=0; k<N; ++k) {
            R   =   fabs(theta(j)-theta(k));
            if (R==0) {
                K(j,k)  =   0;
            }
            else {
                K(j,k)  =   exp(-scale*R*R);
            }
        }
    }

    #elif SINC
    /****************************/
    /*      K   =   sin(R)/R    */
    /****************************/
    cout << "Matrix is sin(R)/R." << endl;
    VectorXd theta  =   VectorXd::Random(M);
    sort(theta.data(), theta.data()+theta.size());
    K   =   MatrixXd(M,N);
    double R;
    double scale    =   1;
    for (unsigned j=0; j<M; ++j) {
        for (unsigned k=0; k<N; ++k) {
            R   =   fabs(theta(j)-theta(k));
            if (R==0) {
                K(j,k)  =   0;
            }
            else {
                K(j,k)  =   sin(scale*R)/R;
            }
        }
    }
    #endif
}

//  FUNCTION:   Obtains an entry of the matrix
double get_Matrix_Entry(const unsigned i, const unsigned j){
    return K(i,j);
}

//  FUNCTION:   Obtains a sub-matrix of the matrix
void get_Matrix(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, MatrixXd& A){
    A   =   K.block(start_Row,start_Col,n_Rows,n_Cols);
}

//  FUNCTION:   Obtains a row of the matrix
void get_Matrix_Row(const unsigned start_Col, const unsigned n_Cols, const unsigned row_Index, VectorXd& v){
    v   =   VectorXd(n_Cols);
    for (unsigned j=0; j<n_Cols; ++j) {
        v(j)    =   get_Matrix_Entry(row_Index,start_Col+j);
    }
}

//  FUNCTION:   Obtains a column of the matrix
void get_Matrix_Col(const unsigned start_Row, const unsigned n_Rows, const unsigned col_Index, VectorXd& v){
    v   =   VectorXd(n_Rows);
    for (unsigned j=0; j<n_Rows; ++j) {
        v(j)    =   get_Matrix_Entry(start_Row+j,col_Index);
    }
}

//  FUNCTION:   Obtains index and value of maximum of the absolute entry in a vector
unsigned max_Abs_Vector(const VectorXd& v, double& max, const vector<int>& not_Allowed_Indices){
    unsigned n      =   v.size();
    max             =   v(0);
    unsigned index  =   0;
    for(unsigned j=0; j<n; ++j){
        if(find(not_Allowed_Indices.begin(),not_Allowed_Indices.end(),j)==not_Allowed_Indices.end()){
            if(fabs(v(j))>fabs(max)){
                max     =   v(j);
                index   =   j;
            }
        }
    }
    return index;
}

//  FUNCTION:   Obtains low-rank decomposition to desired tolerance
void partial_Piv_LU(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const double tolerance, unsigned& computed_Rank, MatrixXd& U, MatrixXd& V){

    //  If the matrix is small enough do not do anything
    unsigned tolerable_Rank =   5;
    if (n_Cols <= tolerable_Rank){
        get_Matrix(start_Row, start_Col, n_Rows, n_Cols,U);
        V               =   MatrixXd::Identity(n_Cols,n_Cols);
        computed_Rank   =   n_Cols;
        return;
    }
    else if (n_Rows <= tolerable_Rank){
        U               =   MatrixXd::Identity(n_Rows,n_Rows);
        get_Matrix(start_Row, start_Col, n_Rows, n_Cols,V);
        computed_Rank   =   n_Rows;
        return;
    }

    vector<int> row_Index;
    vector<int> col_Index;
    vector<VectorXd> u;
    vector<VectorXd> v;

    double max, gamma, unused_max;

    /*  INITIALIZATION  */

    //  Initialize the matrix norm and the the first row index
    double matrix_Norm  =   0;
    row_Index.push_back(0);

    unsigned pivot;

    computed_Rank   =   0;

    VectorXd a, row, col;

    double row_Squared_Norm, row_Norm, col_Squared_Norm, col_Norm;
    //  Repeat till the desired tolerance is obtained
    do {
        //  Generation of the row
        get_Matrix_Row(start_Col, n_Cols, start_Row+row_Index.back(), a);

        //  Row of the residuum and the pivot column
        row =   a;
        for (unsigned l=0; l<computed_Rank; ++l) {
            row =   row-u[l](row_Index.back())*v[l];
        }
        pivot   =   max_Abs_Vector(row, max, col_Index);

        unsigned max_tries  =   20;
        unsigned count      =    0;
        unsigned count1     =    0;

        //  This randomization is needed if in the middle of the algorithm the row happens to be exactly the linear combination of the previous rows.
        while (fabs(max)<1e-16 && count < max_tries) {
            unsigned new_Row_Index;
            row_Index.pop_back();
            do {
                new_Row_Index   =   rand()%n_Rows;
                ++count1;
            } while (find(row_Index.begin(),row_Index.end(),new_Row_Index)!=row_Index.end() && count1 < max_tries);
            count1  =   0;
            row_Index.push_back(new_Row_Index);

            //  Generation of the row
            get_Matrix_Row(start_Col, n_Cols, start_Row+row_Index.back(), a);

            //  Row of the residuum and the pivot column
            row =   a;
            for (unsigned l=0; l<computed_Rank; ++l) {
                row =   row-u[l](row_Index.back())*v[l];
            }
            pivot   =   max_Abs_Vector(row, max, col_Index);
            ++count;
        }
        
        if (count == max_tries) break;

        count = 0;

        col_Index.push_back(pivot);

        //  Normalizing constant
        gamma   =   1.0/max;

        //  Generation of the column
        get_Matrix_Col(start_Row, n_Rows, start_Col+col_Index.back(), a);

        //  Column of the residuum and the pivot row
        col =   a;
        for (unsigned l=0; l<computed_Rank; ++l) {
            col =   col-v[l](col_Index.back())*u[l];
        }
        pivot   =   max_Abs_Vector(col, unused_max, row_Index);


        //  This randomization is needed if in the middle of the algorithm the columns happens to be exactly the linear combination of the previous columns.
        while (fabs(max)<1e-16 && count < max_tries) {
            col_Index.pop_back();
            unsigned new_Col_Index;
            do {
                new_Col_Index   =   rand()%n_Cols;
            } while (find(col_Index.begin(),col_Index.end(),new_Col_Index)!=col_Index.end() && count1 < max_tries);
            count1  =   0;
            col_Index.push_back(new_Col_Index);

            //  Generation of the column
            get_Matrix_Col(start_Row, n_Rows, start_Col+col_Index.back(), a);

            //  Column of the residuum and the pivot row
            col =   a;
            for (unsigned l=0; l<computed_Rank; ++l) {
                col =   col-u[l](col_Index.back())*v[l];
            }
            pivot   =   max_Abs_Vector(col, unused_max, row_Index);
            ++count;
        }

        if (count == max_tries) break;

        count = 0;

        row_Index.push_back(pivot);

        //  New vectors
        u.push_back(gamma*col);
        v.push_back(row);

        //  New approximation of matrix norm
        row_Squared_Norm    =   row.squaredNorm();
        row_Norm            =   sqrt(row_Squared_Norm);

        col_Squared_Norm    =   col.squaredNorm();
        col_Norm            =   sqrt(col_Squared_Norm);

        matrix_Norm         =   matrix_Norm +   gamma*gamma*row_Squared_Norm*col_Squared_Norm;
        
        for (unsigned j=0; j<computed_Rank; ++j) {
            matrix_Norm     =   matrix_Norm +   2.0*(u[j].dot(u.back()))*(v[j].dot(v.back()));
        }
        ++computed_Rank;
    } while (row_Norm*col_Norm > fabs(max)*tolerance*matrix_Norm && computed_Rank <= fmin(n_Rows, n_Cols));

    //  If the computed_Rank is close to full-rank then return the trivial full-rank decomposition
    if (computed_Rank>=fmin(n_Rows, n_Cols)) {
        if (n_Rows < n_Cols) {
            U   =   MatrixXd::Identity(n_Rows,n_Rows);
            get_Matrix(start_Row, start_Col, n_Rows, n_Cols,V);
            computed_Rank   =   n_Rows;
            return;
        }
        else {
            get_Matrix(start_Row, start_Col, n_Rows, n_Cols,U);
            V   =   MatrixXd::Identity(n_Cols,n_Cols);
            computed_Rank   =   n_Cols;
            return;
        }
    }

    U   =   MatrixXd(n_Rows,computed_Rank);
    V   =   MatrixXd(computed_Rank,n_Cols);
    for (unsigned j=0; j<computed_Rank; ++j) {
        U.col(j)    =   u[j];
        V.row(j)    =   v[j];
    }
}

//  FUNCTION:   Performs the fast solve recursively
void fast_Solve(const unsigned start_Row, const unsigned start_Col, const unsigned n_Rows, const unsigned n_Cols, const unsigned n_Leaf, const MatrixXd& b, MatrixXd& x){
    if (n_Rows <= n_Leaf) {
        MatrixXd A;
        get_Matrix(start_Row, start_Col, n_Rows, n_Cols, A);
        x   =   A.fullPivLu().solve(b);
    }
    else{
        unsigned left_Start_Row =   start_Row;
        unsigned left_Rows      =   n_Rows/2;
        unsigned left_Start_Col =   start_Col;
        unsigned left_Cols      =   n_Cols/2;

        unsigned right_Start_Row=   start_Row+left_Rows;
        unsigned right_Rows     =   n_Rows-left_Rows;
        unsigned right_Start_Col=   start_Col+left_Cols;
        unsigned right_Cols     =   n_Cols-left_Cols;

        double low_Rank_Tolerance   =   1e-16;
        unsigned computed_Rank_1, computed_Rank_2, n_Rhs;

        n_Rhs       =   b.cols();

        //  Get the low rank form of the two off diagonal blocks
        MatrixXd U1, V1, U2, V2;

        partial_Piv_LU(left_Start_Row, right_Start_Col, left_Rows, right_Cols, low_Rank_Tolerance, computed_Rank_1, U1, V1);

        partial_Piv_LU(right_Start_Row, left_Start_Col, right_Rows, left_Cols, low_Rank_Tolerance, computed_Rank_2, U2, V2);


        //  Get the right hand side for the two smaller linear systems

        MatrixXd b1(left_Rows,computed_Rank_1+n_Rhs);

        b1.block(0,0,left_Rows,computed_Rank_1)     =   U1;
        b1.block(0,computed_Rank_1,left_Rows,n_Rhs) =   b.block(0,0,left_Rows,n_Rhs);

        MatrixXd b2(right_Rows,computed_Rank_2+n_Rhs);

        b2.block(0,0,right_Rows,computed_Rank_2)    =   U2;
        b2.block(0,computed_Rank_2,right_Rows,n_Rhs)=   b.block(left_Rows,0,right_Rows,n_Rhs);


        //  Obtain the solution for the two smaller linear systems
        MatrixXd Y1, Y2;
        fast_Solve(left_Start_Row, left_Start_Col, left_Rows, left_Cols, n_Leaf, b1, Y1);
        fast_Solve(right_Start_Row, right_Start_Col, right_Rows, right_Cols, n_Leaf, b2, Y2);

        //  Combine the two solutions to obtain the solution of the larger linear system
        MatrixXd Y3 =   V1*Y2;
        MatrixXd Y4 =   V2*Y1;

        
        //  Gets the intermediate smaller matrix which needs to be solver for
        MatrixXd intermediate_Matrix    =   MatrixXd::Identity(computed_Rank_1+computed_Rank_2,computed_Rank_1+computed_Rank_2);

        intermediate_Matrix.block(0, computed_Rank_1, computed_Rank_1, computed_Rank_2) =   Y3.block(0,0,computed_Rank_1,computed_Rank_2);
        intermediate_Matrix.block(computed_Rank_1, 0, computed_Rank_2, computed_Rank_1) =   Y4.block(0,0,computed_Rank_2,computed_Rank_1);

        
        MatrixXd Y(computed_Rank_1+computed_Rank_2,n_Rhs);
        Y.block(0,0,computed_Rank_1,n_Rhs)              =   Y3.block(0,computed_Rank_2,computed_Rank_1,n_Rhs);
        Y.block(computed_Rank_1,0,computed_Rank_2,n_Rhs)=   Y4.block(0,computed_Rank_1,computed_Rank_2,n_Rhs);

        MatrixXd temp_Sol   =   intermediate_Matrix.fullPivLu().solve(Y);

        x   =   MatrixXd(n_Rows, n_Rhs);

        x.block(0,0,left_Rows,n_Rhs)            =   Y1.block(0, computed_Rank_1, left_Rows, n_Rhs) - Y1.block(0, 0, left_Rows, computed_Rank_1)*temp_Sol.block(0, 0, computed_Rank_1, n_Rhs);
        x.block(left_Rows,0,right_Rows,n_Rhs)   =   Y2.block(0, computed_Rank_2, right_Rows, n_Rhs) - Y2.block(0, 0, right_Rows, computed_Rank_2)*temp_Sol.block(computed_Rank_1, 0, computed_Rank_2, n_Rhs);
    }
}

//  FUNCTION: This is the function user will call.
void HODLR_Solve(const MatrixXd& b, const unsigned n_Leaf, MatrixXd& x){
    unsigned N  =   b.rows();
    fast_Solve(0, 0, N, N, n_Leaf, b, x);
}
