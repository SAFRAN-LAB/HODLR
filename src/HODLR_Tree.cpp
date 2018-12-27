#include "HODLR_Tree.hpp"

// Contructor for the HODLR Tree class:
HODLR_Tree::HODLR_Tree(int n_levels, double tolerance, HODLR_Matrix* A) 
{
    this->n_levels  = n_levels;
    this->tolerance = tolerance;
    this->A         = A;
    this->N         = A->N;
    nodes_in_level.push_back(1);

    for (int j = 1; j <= n_levels; j++) 
    {
        nodes_in_level.push_back(2 * nodes_in_level.back());
    }
    
    this->createTree();
}

// Destructor:
HODLR_Tree::~HODLR_Tree() 
{
    for(int j = 0; j <= n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            delete tree[j][k];
        }
    }
}

// Creates the node at the root level:
void HODLR_Tree::createRoot() 
{
    HODLR_Node* root = new HODLR_Node(0, 0, 0, 0, N, tolerance);
    std::vector<HODLR_Node*> level;
    level.push_back(root);
    tree.push_back(level);
}

// Function that adds on children for the given level and node number:
void HODLR_Tree::createChildren(int level_number, int node_number) 
{
    //  Adding left child:
    HODLR_Node* left = new HODLR_Node(level_number + 1, 2 * node_number, 0, 
                                      tree[level_number][node_number]->c_start[0], 
                                      tree[level_number][node_number]->c_size[0], 
                                      tolerance
                                     );
    tree[level_number + 1].push_back(left);

    //  Adding right child
    HODLR_Node* right = new HODLR_Node(level_number + 1, 2 * node_number + 1, 1, 
                                       tree[level_number][node_number]->c_start[1], 
                                       tree[level_number][node_number]->c_size[1], 
                                       tolerance
                                      );
    tree[level_number + 1].push_back(right);
}

// Creates the tree structure:
// Depending upon the parameters set by the user, this function
// creates a tree having required number of nodes in the tree
void HODLR_Tree::createTree() 
{
    this->createRoot();
    
    for(int j = 0; j < n_levels; j++) 
    {
        std::vector<HODLR_Node*> level;
        tree.push_back(level);
        
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            this->createChildren(j, k);
        }
    }
}

// This function assembles the elements of the tree:
// That is:
// For leaf nodes, it directly evaluates the matrix entries:
// For nonleaf nodes, it gets the low rank representation of the underlying matrix:
void HODLR_Tree::assembleTree(bool is_sym, bool is_pd) 
{
    this->is_sym = is_sym;
    this->is_pd  = is_pd;
    // Assembly of nonleaf nodes:
    for(int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            tree[j][k]->assembleNonLeafNode(A, is_sym);
        }
    }

    // Assembly of leaf nodes:
    #pragma omp parallel for
    for (int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        tree[n_levels][k]->assembleLeafNode(A);
    }
}

// Used mainly to debug:
// Prints the details of the considered nodes:
void HODLR_Tree::printNodeDetails(int level_number, int node_number)
{
    tree[level_number][node_number]->printNodeDetails();
}

// Used mainly to debug:
// Prints details of all the nodes in the tree:
void HODLR_Tree::printTreeDetails()
{
    for(int j = 0; j <= n_levels; j++)
    {
        for(int k = 0; k < nodes_in_level[j]; k++)
        {
            this->printNodeDetails(j, k);
            cout << "=======================================================================================================================================" << endl;
        }
        cout << endl << endl;
    }
}

// Performs a MatMat product with the given matrix X.
Mat HODLR_Tree::matmatProduct(Mat x) 
{
    // Initializing matrix b:
    Mat b = Mat::Zero(N, x.cols());

    // At non-leaf levels:
    for (int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            tree[j][k]->matmatProductNonLeaf(x, b);
        }
    }

    // At leaf level:
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        tree[n_levels][k]->matmatProductLeaf(x, b);
    }

    return b;
} 

// Factorizes the matrix to get LU for each block in case of the nonSPD
// and L for each block in case the matrix is SPD
void HODLR_Tree::factorize()
{
    if(is_sym == true && is_pd == true)
        this->factorizeSPD();
    else
        this->factorizeNonSPD();
}

// Returns x, by solving Ax = b, where A is the matrix represented by the HODLR structure
Mat HODLR_Tree::solve(Mat b)
{
    if(is_sym == true && is_pd == true)
        return this->solveSPD(b);
    else
        return this->solveNonSPD(b);
}

// Returns the log determinant of the matrix represented by the HODLR structure
dtype HODLR_Tree::logDeterminant()
{
    if(is_sym == true && is_pd == true)
        return this->logDeterminantSPD();
    else
        return this->logDeterminantNonSPD();
}

// Creates a plot of the HODLR matrix represented.
// Basically, this function shows the extent of "lowrankness" through
// different intensity of colors in the plots.
// Additionally, it also shows the ranks of the blocks. Again useful to debug:
void HODLR_Tree::plotTree()
{
    std::string HODLR_PATH = std::getenv("HODLR_PATH");
    std::string FILE       = HODLR_PATH + "/src/plot_tree.py";
    std::string COMMAND    = "python " + FILE;
    std::ofstream myfile;
    myfile.open ("rank.txt");
    // First entry is the number of levels:
    myfile << n_levels << endl;
    for(int j = 0; j < n_levels; j++)
    {
        for(int k = 0; k < nodes_in_level[j]; k++)
        {
            myfile << tree[j][k]->rank[0] << endl;
        }
    }
    // Closing the file:
    myfile.close();
    cout << "Plotting Tree..." << endl;
    system(COMMAND.c_str());
}
