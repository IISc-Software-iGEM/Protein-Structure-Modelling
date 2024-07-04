#include <bits/stdc++.h>
#include "../../Mosek/matrixoperations.h"
#include "../../../eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;

void saveMatrix(const MatrixXd& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << matrix.rows() << " " << matrix.cols() << "\n";
        file << matrix << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

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

vector<MatrixXd> extractRotations(const MatrixXd& W, int dim, int num_patch) {
    int Md = W.cols();  // Assuming W is of size Md x Md
    // Initialize ortho_rot and rot
    MatrixXd ortho_rot = MatrixXd::Zero(dim, Md);
    std::vector<MatrixXd> rot(num_patch);

    // Loop through patches
    for (int i = 0; i < num_patch; ++i) {
        // Select the i-th dim columns of W
        MatrixXd W_patch = W.block(0, dim * i, W.rows(), dim);
        
        // Perform SVD on the selected columns
        JacobiSVD<MatrixXd> svd(W_patch, ComputeThinU | ComputeThinV);
        MatrixXd U = svd.matrixU();
        MatrixXd V = svd.matrixV();
        
        // Compute U * V'
        MatrixXd UVT = U * V.transpose();
        
        // Store in ortho_rot
        ortho_rot.block(0, dim * i, dim, dim) = UVT;
        
        // Store in rot
        rot[i] = UVT;
    }

    return rot;
}

int main() {
    string filename = "matrix.txt";
    
    // Load the matrix from the file
    MatrixXd B = loadMatrix(filename);

    int row = B.rows();
    int col = B.cols();

    int dim = 3;
    int Md = 24;
    int num_patch = 8; // Number of patches, assuming num_patch * dim = Md
    cout << row << " " << col << " " << endl;

    // Eigenvalue decomposition (from previous code)
    EigenSolver<MatrixXd> eigensolver(B);
    VectorXd eigenvalues = eigensolver.eigenvalues().real();
    MatrixXd eigenvectors = eigensolver.eigenvectors().real();

    std::vector<std::pair<double, VectorXd>> eigenPairs;
    for (int i = 0; i < Md; i++) {
        eigenPairs.push_back(std::make_pair(eigenvalues(i), eigenvectors.col(i)));
    }
    std::sort(eigenPairs.begin(), eigenPairs.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first > rhs.first; // Sort in descending order
    });

    MatrixXd V(Md, dim);
    VectorXd S(dim);
    for (int i = 0; i < dim; i++) {
        S(i) = eigenPairs[i].first;
        V.col(i) = eigenPairs[i].second;
    }

    MatrixXd S_diag = S.asDiagonal();
    MatrixXd sqrt_S_diag = S_diag.cwiseSqrt();
    MatrixXd W = sqrt_S_diag * V.transpose();

    cout << W << endl;

    // Extract rotations
    vector<MatrixXd> rotations = extractRotations(W, dim, num_patch);

    // Output results
    for (int i = 0; i < rotations.size(); ++i) {
        cout << "Rotation matrix " << i + 1 << ":\n" << rotations[i] << endl;
    }

    // Convert rotations to a single 3*24 matrix
    MatrixXd R(dim, Md);
    for (int i = 0; i < num_patch; ++i) {
        R.block(0, dim * i, dim, dim) = rotations[i];
    }
    cout << R << endl;
    saveMatrix(R, "rotations.txt");
    return 0;
}
