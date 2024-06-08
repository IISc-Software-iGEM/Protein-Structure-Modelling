#include <bits/stdc++.h>
using namespace std;

vector<vector<double>> transpose(vector<vector<double>> matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    vector<vector<double>> Tmatrix(cols, vector<double>(rows,0));

    for(int i = 0;i < rows; i++) {
        for(int j = 0; j < cols; j++) {
            Tmatrix[j][i] = matrix[i][j];
        }
    }

    return Tmatrix;
}

vector<vector<double>> matrixMult(vector<vector<double>> matrixA, vector<vector<double>> matrixB) {
    int rowsMatrixA = matrixA.size();
    int colsMatrixA = matrixA[0].size();

    int rowsMatrixB = matrixB.size();
    int colsMatrixB = matrixB[0].size();

    // vector to return
    vector< vector <double>> newMatrix(rowsMatrixA, vector<double> (colsMatrixB,0));

    // fill entries
    for(int i = 0;i < rowsMatrixA; i++) {
        for(int j = 0;j < colsMatrixB; j++) {
            for(int k = 0;k < colsMatrixA; k++) {
                newMatrix[i][j] += matrixA[i][k]*matrixB[k][j];
            }
        }
    } 

    return newMatrix;
}