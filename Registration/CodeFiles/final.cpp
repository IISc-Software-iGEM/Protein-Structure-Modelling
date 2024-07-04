#include <bits/stdc++.h>
#include "../../Mosek/matrixoperations.h"
#include "../../../eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;


MatrixXd loadMatrix(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Unable to open file for reading: " << filename << std::endl;
        return MatrixXd();
    }

    int rows, cols;
    file >> rows >> cols;

    MatrixXd matrix(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file >> matrix(i, j);
        }
    }
    
    file.close();
    return matrix;
}

int main() {
    // Define file paths and dimensions
    string filename_B = "B_eigen.txt";
    string filename_L = "L_eigen.txt";
    string filename_ortho_rot = "rotations.txt";

    int dim = 3;        // Dimension of each block in ortho_rot
    int Md = 24;        // Number of columns in B and rows in L
    int num_pts = 1031; // Number of points to select from x_cord

    // Load matrices B, L, and ortho_rot from files
    MatrixXd B = loadMatrix(filename_B);
    MatrixXd L = loadMatrix(filename_L);
    MatrixXd ortho_rot = loadMatrix(filename_ortho_rot);

    // Compute L_inv = pinv(L)
    MatrixXd L_inv = L.completeOrthogonalDecomposition().pseudoInverse();

    for(int i = 0;i < ortho_rot.cols();i++) {
        cout << ortho_rot(0,i) << " ";
    }
    cout << endl;

    // Compute x_cord = ortho_rot * B * L_inv
    MatrixXd x_cord = ortho_rot * B * L_inv;

    // Select first num_pts columns of x_cord
    x_cord.conservativeResize(x_cord.rows(), num_pts);

    // Output results column-wise
    for (int i = 0; i < x_cord.cols(); ++i) {
        cout << x_cord.col(i).transpose() << endl;
    }

    return 0;
}