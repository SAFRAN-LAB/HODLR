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

void HODLR_Tree::createRoot() 
{
    HODLR_Node* root = new HODLR_Node(0, 0, 0, 0, N, tolerance);
    std::vector<HODLR_Node*> level;
    level.push_back(root);
    tree.push_back(level);
}

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

void HODLR_Tree::assembleTree(bool is_spd, bool is_sym) 
{
    this->is_spd = is_spd;
    this->is_sym = is_sym;
    // Assembly of nonleaf nodes:
    for(int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for(int k = 0; k < nodes_in_level[j]; k++) 
        {
            tree[j][k]->assembleNonLeafNode(A, is_spd, is_sym);
        }
    }

    // Assembly of leaf nodes:
    #pragma omp parallel for
    for (int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        tree[n_levels][k]->assembleLeafNode(A);
    }
}

void HODLR_Tree::printNodeDetails(int level_number, int node_number)
{
    tree[level_number][node_number]->printNodeDetails();
}

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

void HODLR_Tree::matmatProduct(MatrixXd x, MatrixXd& b) 
{
    // Initializing matrix b:
    b = MatrixXd::Zero(N, x.cols());

    // At non-leaf levels:
    for (int j = 0; j < n_levels; j++) 
    {
        #pragma omp parallel for
        for (int k = 0; k < nodes_in_level[j]; k++) 
        {
            tree[j][k]->matmatProductNonLeaf(x, b, is_spd);
        }
    }

    // At leaf level:
    #pragma omp parallel for
    for(int k = 0; k < nodes_in_level[n_levels]; k++) 
    {
        tree[n_levels][k]->matmatProductLeaf(x, b);
    }
} 

void HODLR_Tree::factorize()
{
    if(is_spd == true)
        this->factorizeSPD();
    else
        this->factorizeNonSPD();
}

MatrixXd HODLR_Tree::solve(MatrixXd b)
{
    if(is_spd == true)
        return this->solveSPD(b);
    else
        return this->solveNonSPD(b);
}

double HODLR_Tree::logDeterminant()
{
    if(is_spd == true)
        return this->logDeterminantSPD();
    else
        return this->logDeterminantNonSPD();
}

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
