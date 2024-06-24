#include <bits/stdc++.h>
#include "../boundreader.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <chrono>
// #include "fusion.h"

using namespace std;
// using namespace Eigen;
// using namespace mosek::fusion;
// using namespace monty;

int main() {

    auto start = std::chrono::high_resolution_clock::now();
    int dimension = 41;

    string equalityBoundFile = "eqbnds.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityBoundFile);

    string upperBoundFile = "ubnds.txt";
    vector < boundReader* > UpperBounds = bound_reader(upperBoundFile);

    int n = EqualityBounds.size();
    int m = UpperBounds.size();

    cout << "Number of equality bounds: " << n << endl;
    cout << "Number of upper bounds: " << m << endl;

    // reading U matrix
    int k = 17;
    Eigen::MatrixXd U(dimension, k);

    ifstream inputfile("U.txt");
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
    
    cout << "U dimensions: ";
    cout << "rows: " << U.rows() << " cols: " << U.cols() << endl;

    // Sum the columns of U to get a row vector
    Eigen::RowVectorXd colSum = U.colwise().sum();
    cout << "ColSum: " << colSum << endl;

    // Transpose to get a column vector
    Eigen::VectorXd UTe = colSum.transpose();
    cout << "UTe: " << UTe << endl;

    // Perform QR decomposition
    Eigen::HouseholderQR<Eigen::MatrixXd> qr(UTe);
    Eigen::MatrixXd V = qr.householderQ();

    // Update U by multiplying with V(:,2:end)
    Eigen::MatrixXd UV = U * V.block(0, 1, V.rows(), V.cols() - 1);

    cout << "dimension of UV=> ";
    cout << "rows: " << UV.rows() << " cols: " << UV.cols() << endl;
    // data for vec cons

    // Model::t M = new Model("ayushRaina");
    // auto _M = finally([&]() {M->dispose();});

    // Variable::t XX = M->variable(Domain::inPSDCone(k-1));
    // auto traceX = Expr::sum(XX->diag());

    cout << "Data for vector constraints:" << endl;
    for(int i = 0;i < EqualityBounds.size(); i++) {
        boundReader* X = EqualityBounds[i];
        Eigen::MatrixXd A = UV.row(X->x_i-1) - UV.row(X->x_j-1);

        auto XX = A.transpose() * A;
        cout << XX << endl;

        auto reqVector = new_array_ptr<double,2> (shape(k-1, k-1));
        for(int i = 0;i < k-1; i++) {
            for(int j = 0;j < k-1; j++) {
                reqVector->operator()(i,j) = XX(i,j);
            }
        } 

        break;
    }

    // cout << "Data for b: " << endl;
    // for(int i = 0;i < EqualityBounds.size(); i++) {
    //     boundReader* X = EqualityBounds[i];
    //     cout << pow(X->bound,2) << " ";
    // }
    cout << endl;

    cout << "Data Part II: " << endl;

    // for(int i = 0;i < UpperBounds.size(); i++) {
    //     boundReader* X = UpperBounds[i];
    //     cout << UV.row(X->x_i-1) - UV.row(X->x_j-1) << endl;
    // }

    // cout << "Data for b: " << endl;
    // for(int i = 0;i < UpperBounds.size(); i++) {
    //     boundReader* X = UpperBounds[i];
    //     cout << pow(X->bound,2) << " ";
    // }

    cout << "dimensions of V=> ";
    cout << "rows: " << V.rows() << " cols: " << V.cols() << endl;


    // Record end time in milliseconds
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Time taken: " << duration.count() << " milliseconds" << endl;

}