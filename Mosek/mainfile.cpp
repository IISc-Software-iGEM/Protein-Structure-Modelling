#include "matrixoperations.h"
#include "boundreader.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;
using namespace std;
using namespace Eigen;

int calculateRank(const MatrixXd& matrix) {
    FullPivLU< MatrixXd > lu_decomp(matrix);
    return lu_decomp.rank();
}

Eigen::MatrixXd approximateRank(const Eigen::MatrixXd& matrix, int rank) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();

    // Truncate to the desired rank
    Eigen::MatrixXd truncatedU = U.leftCols(rank);
    Eigen::MatrixXd truncatedV = V.leftCols(rank);
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(rank, rank);

    for (int i = 0; i < rank; ++i) {
        S(i, i) = singularValues(i);
    }

    return truncatedU * S * truncatedV.transpose();
}

// Function to extract 3D coordinates using Cholesky decomposition
Eigen::MatrixXd extract3DCoordinates(const Eigen::MatrixXd& matrix) {
    // Perform Cholesky decomposition
    Eigen::LLT<Eigen::MatrixXd> llt(matrix);
    Eigen::MatrixXd L = llt.matrixL();
    
    // Extract the first 3 columns for 3D coordinates
    return L.leftCols(3);
}

int main() {
    int dimension = 133;
    int numConstraints = 1;
    
    string equalityBoundFile = "testfile.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityBoundFile);

    int n = EqualityBounds.size();
    cout << "Total number of equality bounds are: " << n << endl; 

    Model::t M = new Model("ayushRaina");
    auto _M = finally([&]() {M->dispose();});

    Variable::t XX = M->variable(Domain::inPSDCone(dimension));
    auto traceX = Expr::sum(XX->diag());

    //  Set Objective Function
    M->objective(ObjectiveSense::Minimize, traceX);
    
    // setting equality bounds
    for(int i = 0;i < n; i++) {
        boundReader* X = EqualityBounds[i];

        // create E
        vector<vector< double >> E(1, vector<double>(dimension,0));
        E[0][X->x_i-1] = 1;
        E[0][X->x_j-1] = -1;

        // E transpose
        vector<vector <double> > ETranspose = transpose(E);

        // matmul
        vector<vector<double>> answer = matrixMult(ETranspose, E);

        int row = answer.size();
        int col = answer[0].size();

        shared_ptr<ndarray<double, 2>> A = new_array_ptr<double, 2>(shape(row,col));
        for(int i = 0;i < n; i++) { 
            for(int j = 0;j < n; j++) {
                A->operator()(i, j) = answer[i][j];
            }
        }

        auto EZ = Expr::mulDiag(A,XX);
        auto traceEZ = Expr::sum(EZ);

        M->constraint(traceEZ, Domain::equalsTo(X->bound));
    }

    // Setting upper bounds
    string upperBoundsFile = "upBounds.txt";
    vector < boundReader* > UpperBounds = bound_reader(upperBoundsFile);
    int m = Y.size();

    cout << "Total number of upper bounds are: " << m << endl;

    for(int i = 0;i < m; i++) {
        boundReader* Y = UpperBounds[i];

        // create E
        vector<vector< double >> E(1, vector<double>(dimension,0));
        E[0][Y->x_i-1] = 1;
        E[0][Y->x_j-1] = -1;

        // E transpose
        vector<vector <double> > ETranspose = transpose(E);

        // matmul
        vector<vector<double>> answer = matrixMult(ETranspose, E);

        int row = answer.size();
        int col = answer[0].size();

        shared_ptr<ndarray<double, 2>> A = new_array_ptr<double, 2>(shape(row,col));
        for(int i = 0;i < row; i++) { 
            for(int j = 0;j < col; j++) {
                A->operator()(i, j) = answer[i][j];
            }
        }

        auto EZ = Expr::mulDiag(A,XX);
        auto traceEZ = Expr::sum(EZ);

        M->constraint(traceEZ, Domain::lessThan(X->bound));
    }

    M->solve();
    cout << "Solution : " << endl;
    
    // making eigen matrix;
    auto XXsolution = XX->level();
    MatrixXd SolutionMatrix(dimension, dimension);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            SolutionMatrix(i, j) = XXsol->operator[](i * dimension + j);
        }
    }

    // Print the solution matrix
    cout << "Solution = \n" << SolutionMatrix << endl;

    // Calculate the rank of the solution matrix
    int rank = calculateRank(SolutionMatrix);
    cout << "Rank of the solution matrix: " << rank << endl;

    // Approximate the solution matrix to rank-3
    Eigen::MatrixXd rank3Matrix = approximateRank(solutionMatrix, 3);
    cout << "Rank-3 approximation of the solution matrix: \n" << rank3Matrix << endl;

    // Extract 3D coordinates using Cholesky decomposition
    Eigen::MatrixXd coordinates = extract3DCoordinates(rank3Matrix);
    cout << "3D coordinates of the atoms: \n" << coordinates << endl;

    return 0;
}