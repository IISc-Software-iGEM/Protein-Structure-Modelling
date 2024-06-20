#include <bits/stdc++.h>
#include "../boundreader.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/Sparse"
#include <chrono>

using namespace std;

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
    
    


    // Record end time in milliseconds
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Time taken: " << duration.count() << " milliseconds" << endl;

}