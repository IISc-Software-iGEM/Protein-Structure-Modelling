#include "matrixoperations.h"
#include "boundreader.h"
#include "fusion.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
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

Eigen::MatrixXd getUVmatrix(const string& filename, int k, int dimension) {
    Eigen::MatrixXd U(dimension, k);

    ifstream inputfile(filename);
    if(!inputfile.is_open()) {
        cout << "Error opening file" << endl;
        exit(0);
    }

    for(int i = 0; i < dimension; i++) {
        for(int j = 0; j < k; j++) {
            inputfile >> U(i,j);
        }
    }
    inputfile.close();

    // Sum the columns of U to get a row vector
    Eigen::RowVectorXd colSum = U.colwise().sum();

    // Transpose to get a column vector
    Eigen::VectorXd UTe = colSum.transpose();

    // Perform QR decomposition
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(UTe);
    Eigen::MatrixXd V = qr.householderQ();

    // Update U by multiplying with V(:,2:end)
    Eigen::MatrixXd UV = U * V.block(0, 1, V.rows(), V.cols() - 1);

    cout << "dimension of UV=> ";
    cout << "rows: " << UV.rows() << " cols: " << UV.cols() << endl;

    cout << "dimensions of V=> ";
    cout << "rows: " << V.rows() << " cols: " << V.cols() << endl;

    return UV;
}

int main() {

    // Record start time
    auto start = std::chrono::high_resolution_clock::now();

    int dimension = 133;
    int k = 61;
    string equalityBoundFile = "eqbnds.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityBoundFile);

    int n = EqualityBounds.size();
    cout << "Total number of equality bounds are: " << n << endl;  

    // get UV matrix
    Eigen::MatrixXd UV = getUVmatrix("U.txt", k, dimension);

    Model::t M = new Model("ayushRaina");
    auto _M = finally([&]() {M->dispose();});

    Variable::t XX = M->variable(Domain::inPSDCone(k-1));
    auto traceX = Expr::sum(XX->diag());
    
    // setting equality bounds
    for(int i = 1;i < n; i+=4) {
        boundReader* X = EqualityBounds[i];

        // creating A_ij dash
        // creating A_ij dash
        Eigen::MatrixXd A_new = UV.row(X->x_i - 1) - UV.row(X->x_j - 1);
        auto A_ijDash = A_new.transpose() * A_new;
        cout << A_ijDash.rows() << " " << A_ijDash.cols() << endl;

        // create vector of A_ij dash
        auto reqVector = new_array_ptr<double,2> (shape(k-1, k-1));
        for(int i = 0;i < k-1; i++) {
            for(int j = 0;j < k-1; j++) {
                reqVector->operator()(i,j) = A_ijDash(i,j);
            }
        }

        auto reqMatrix = mosek::fusion::Matrix::dense(reqVector);
        auto Imatrix = mosek::fusion::Matrix::eye(k-1);

        auto traceEZ = Expr::dot(reqMatrix,XX);
	    cout << "matmul done for " << i << "th time" << endl;
        M->constraint(traceEZ, Domain::equalsTo(X->bound));
	    cout << i << "th equality constraint added" << endl;
    }

    // Setting upper bounds
    string upperBoundsFile = "ubnds.txt";
    vector < boundReader* > UpperBounds = bound_reader(upperBoundsFile);
    int m = UpperBounds.size();
    Variable::t xi_1 = M->variable(m, Domain::greaterThan(0.0));
    cout << "Total number of upper bounds are: " << m << endl;

    for(int i = 0;i < m; i++) {
        boundReader* Y = UpperBounds[i];

        // creating A_ij dash
        Eigen::MatrixXd A_new = UV.row(Y->x_i - 1) - UV.row(Y->x_j - 1);

        auto A_ijDash = A_new.transpose() * A_new;
        cout << A_ijDash.rows() << " " << A_ijDash.cols() << endl;

        // create vector of A_ij dash
        auto reqVector = new_array_ptr<double,2> (shape(k-1, k-1));
        for(int i = 0;i < k-1; i++) {
            for(int j = 0;j < k-1; j++) {
                reqVector->operator()(i,j) = A_ijDash(i,j);
            }
        }
        
        auto reqMatrix = mosek::fusion::Matrix::dense(reqVector);
        auto Imatrix = mosek::fusion::Matrix::eye(k-1);

        auto traceEZ = Expr::dot(reqMatrix,XX);
        cout << "matmul done for the " << i << "th time" << endl;

        M->constraint(Expr::sub(traceEZ,xi_1->index(i)), Domain::lessThan(Y->bound));
	    cout << i << "th upperBound constraint added" << endl;
    }

    int wIJ = 300;
    int wIJdash = 300;

    auto o2 = Expr::mul(1,xi_1->index(0));
    for(int i = 1;i < m; i++) {
        o2 = Expr::add(o2,Expr::mul(205,xi_1->index(i)));
    }
    //  Set Objective Function
    //M->objective(ObjectiveSense::Minimize, Expr::add(Expr::mul(10,traceX), Expr::add(Expr::mul(wIJ,Expr::sum(xi_1)), Expr::mul(wIJdash,Expr::sum(xi_2)))));
    M->objective(ObjectiveSense::Minimize, Expr::add(Expr::mul(1,traceX), o2));
    cout << "Objective function set" << endl;
    M->setLogHandler([ = ](const std::string & msg) { std::cout << msg << std::flush;} );
    M->solve();
    cout << "Solution : ayush" << endl;
    
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

    //Approximate the solution matrix to rank-3
    Eigen::MatrixXd rank3Matrix = approximateRank(SolutionMatrix, 3);
    //cout << "Rank-3 approximation of the solution matrix: \n" << rank3Matrix << endl;

    rank = calculateRank(rank3Matrix);
    cout << "Rank of the solution matrix: " << rank << endl;

    //Extract 3D coordinates using Cholesky decomposition
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

