#include "matrixoperations.h"
#include "boundreader.h"
#include "fusion.h"
#include "eigen-3.4.0/Eigen/Dense"
#include <chrono>

using namespace mosek::fusion;
using namespace monty;
using namespace std;
using namespace Eigen;

int calculateRank(const MatrixXd& matrix) {
    FullPivLU< MatrixXd > lu_decomp(matrix);
    return lu_decomp.rank();
}

Eigen::MatrixXd approximateRank(const Eigen::MatrixXd& A, int rank) {
    // Eigenvalue decomposition
    SelfAdjointEigenSolver<MatrixXd> es(A);
    if (es.info() != Success) {
        cout << "Eigenvalue decomposition failed!" << endl;
        exit(0);
    }

    // Get eigenvalues and eigenvectors
    VectorXd eigenvalues = es.eigenvalues();
    MatrixXd eigenvectors = es.eigenvectors();

    // Retain top 3 eigenvalues and corresponding eigenvectors
    VectorXd top_eigenvalues = eigenvalues.tail(rank);
    MatrixXd top_eigenvectors = eigenvectors.rightCols(rank);

    // Reconstruct the matrix
    MatrixXd reconstructed_matrix = top_eigenvectors * top_eigenvalues.asDiagonal() * top_eigenvectors.transpose();
    return reconstructed_matrix;
}

// Function to extract 3D coordinates using Cholesky decomposition
Eigen::MatrixXd extract3DCoordinates(const Eigen::MatrixXd& reconstructed_matrix) {
    int rank = 3;
    JacobiSVD<MatrixXd> svd(reconstructed_matrix, ComputeThinU | ComputeThinV);
    MatrixXd U = svd.matrixU();
    VectorXd S = svd.singularValues();

    // Extract the 3D coordinates
    MatrixXd coordinates = U.leftCols(rank) * S.head(rank).asDiagonal();
    return coordinates;
}

int main() {

    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    int dimension = 1547;
    
    string equalityBoundFile = "2m4keqbounds.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityBoundFile);

    int n = EqualityBounds.size();
    cout << "Total number of equality bounds are: " << n << endl; 

    Model::t M = new Model("ayushRaina");
    auto _M = finally([&]() {M->dispose();});

    Variable::t XX = M->variable(Domain::inPSDCone(dimension));
    auto traceX = Expr::sum(XX->diag());
    
    // setting equality bounds
    for(int i = 0;i < n; i++) {
        boundReader* X = EqualityBounds[i];

        auto rows = new_array_ptr<int, 1>( {int((X->x_i*1000)-1), int((X->x_j*1000)-1), int((X->x_i*1000)-1), int((X->x_j*1000)-1)} );
	auto cols = new_array_ptr<int, 1>( {int((X->x_i*1000)-1), int((X->x_j*1000)-1), int((X->x_j*1000)-1), int((X->x_i*1000)-1)} );
	
	auto values = new_array_ptr<double, 1>( {1.0,1.0,-1.0,-1.0} );
	auto A = mosek::fusion::Matrix::sparse(dimension, dimension, rows, cols, values);
        
	
        auto traceEZ = Expr::sum(Expr::mulDiag(A,XX));
	cout << "matmul done for " << i << "th time" << endl;
        M->constraint(traceEZ, Domain::equalsTo(X->bound*1000));
	cout << i << "th equality constraint added" << endl;
    }

    Setting upper bounds
    string upperBoundsFile = "2m4kupbounds.txt";
    vector < boundReader* > UpperBounds = bound_reader(upperBoundsFile);
    int m = UpperBounds.size();
    Variable::t xi_1 = M->variable(m, Domain::greaterThan(0.0));
    cout << "Total number of upper bounds are: " << m << endl;

    for(int i = 0;i < m; i++) {
        boundReader* Y = UpperBounds[i];
        
        auto rows = new_array_ptr<int, 1>( {int((Y->x_i*1000)-1), int((Y->x_j*1000)-1), int((Y->x_i*1000)-1), int((Y->x_j*1000)-1)} );
        auto cols = new_array_ptr<int, 1>( {int((Y->x_i*1000)-1), int((Y->x_j*1000)-1), int((Y->x_j*1000)-1), int((Y->x_i*1000)-1)} );

        auto values = new_array_ptr<double, 1>( {1.0,1.0,-1.0,-1.0} );
        auto A = mosek::fusion::Matrix::sparse(dimension, dimension, rows, cols, values);
        
	auto traceEZ = Expr::sum(Expr::mulDiag(A,XX));
	cout << "matmul done for the " << i << "th time" << endl;

    M->constraint(Expr::sub(traceEZ,xi_1->index(i)), Domain::lessThan(Y->bound*1000));
	cout << i << "th upperBound constraint added" << endl;
    }

    string lowerBoundFile = "2m4klobounds.txt";
    vector < boundReader* > LowerBounds = bound_reader(lowerBoundFile);
    int o = LowerBounds.size();
    cout << "Total Number of lower bounds are: " << o << endl;
    Variable::t xi_2 = M->variable(o, Domain::greaterThan(0.0));
    
    for(int i = 0;i < o; i++) {
       boundReader* Z = LowerBounds[i];

        auto rows = new_array_ptr<int, 1>( {int((Z->x_i*1000)-1), int((Z->x_j*1000)-1), int((Z->x_i*1000)-1), int((Z->x_j*1000)-1)} );
        auto cols = new_array_ptr<int, 1>( {int((Z->x_i*1000)-1), int((Z->x_j*1000)-1), int((Z->x_j*1000)-1), int((Z->x_i*1000)-1)} );

        auto values = new_array_ptr<double, 1>( {1.0,1.0,-1.0,-1.0} );
        auto A = mosek::fusion::Matrix::sparse(dimension, dimension, rows, cols, values);
        
        auto traceEZ = Expr::sum(Expr::mulDiag(A,XX));
        M->constraint(Expr::add(traceEZ,xi_2->index(i)), Domain::greaterThan(Z->bound*1000));
	   cout << i << "th lower bound constraint added" << endl;
    }


    int wIJ = 300;
    int wIJdash = 300;
    //  Set Objective Function
    
    M->objective(ObjectiveSense::Minimize, Expr::add(Expr::mul(0,traceX), Expr::add(Expr::mul(wIJ,Expr::sum(xi_1)), Expr::mul(wIJdash,Expr::sum(xi_2)))));
    
    //M->objective(ObjectiveSense::Minimize, Expr::add(Expr::mul(10,traceX),Expr::mul(wIJ,Expr::sum(xi_1))));

    //M->objective(ObjectiveSense::Minimize, Expr::mul(10,traceX));

    M->setLogHandler([ = ](const std::string & msg) { std::cout << msg << std::flush;} );
    M->solve();
    cout << "Solution : " << endl;
    
    // making eigen matrix;
    auto XXsolution = XX->level();
    MatrixXd SolutionMatrix(dimension, dimension);

    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < dimension; j++) {
            SolutionMatrix(i, j) = XXsolution->operator[](i * dimension + j);
        }
    }

    // Print the solution matrix
    cout << "Solution = \n" << SolutionMatrix << endl;

    // Calculate the rank of the solution matrix
    int rank = calculateRank(SolutionMatrix);
    cout << "Rank of the solution matrix: " << rank << endl;

    // Approximate the solution matrix to rank-3
     Eigen::MatrixXd rank3Matrix = approximateRank(SolutionMatrix, 3);
    //cout << "Rank-3 approximation of the solution matrix: \n" << rank3Matrix << endl;

    rank = calculateRank(rank3Matrix);
    cout << "Rank of the solution matrix: " << rank << endl;

    // Extract 3D coordinates using Cholesky decomposition
    Eigen::MatrixXd coordinates = extract3DCoordinates(rank3Matrix);
    cout << "3D coordinates of the atoms: \n" << coordinates << endl;
    
    // Record end time
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    // Output the duration
    std::cout << "Time taken: " << duration.count() << " milliseconds" << std::endl;
    return 0;
}

