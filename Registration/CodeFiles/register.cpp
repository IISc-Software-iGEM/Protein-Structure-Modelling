#include <bits/stdc++.h>
#include "matrixoperations.h"
#include "eigen-3.4.0/Eigen/Dense"
#include "fusion.h"

using namespace std;
using namespace Eigen;
using namespace monty;
using namespace mosek::fusion;

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

int main() {
    
    auto startTime = chrono::high_resolution_clock::now();
    string filename = "objectiveMatrix.txt";
    
    // Load the matrix from the file
    MatrixXd B = loadMatrix(filename);

    int row = B.rows();
    int col = B.cols();

    shared_ptr<ndarray<double, 2>> A = new_array_ptr<double, 2>(shape(row,col));
    for(int i = 0;i < row; i++) { 
        for(int j = 0;j < col; j++) {
            A->operator()(i, j) = B(i,j);
        }
    }

    Model::t M = new Model("ayushraina");
    auto _M = finally([&]() { M->dispose(); });

    // Creating Semi Definite Variable
    Variable::t X = M->variable("X", Domain::inPSDCone(B.rows()));


    auto AX = Expr::mulDiag(A, X);
    auto trace = Expr::sum(AX);

    // Objective
    M->objective(ObjectiveSense::Maximize, trace); 
    
    int numPatches = 8;
    int dimension = 3;

    for(int i = 0;i < numPatches; i++) {
        auto startIdx = i * dimension;
        auto endIdx = (i + 1) * dimension;

        cout << startIdx << " " << endIdx << endl;
        M->constraint(X->slice(new_array_ptr<int,1>({startIdx,startIdx}), new_array_ptr<int,1>({endIdx,endIdx})), Domain::equalsTo(mosek::fusion::Matrix::eye(dimension)));
    }

    M->setSolverParam("mioTolRelGap", 1e-4);
    M->setSolverParam("mioTolAbsGap", 0.0);

    M->solve();
    std::cout << "Solution : " << std::endl;


    auto Solution = X->level();
    for(int i = 0;i < 24; i++) {
        for(int j = 0;j < 24; j++) {
            B(i,j) = (*Solution)[24*i + j];
        }
    }
    
    saveMatrix(B,"output.txt");
    auto endTime = chrono::high_resolution_clock::now(); 
    cout << "Time Taken: " << chrono::duration_cast<chrono::microseconds>(endTime-startTime).count() << " microseconds" << endl;
    cout << *Solution << endl;  
}
