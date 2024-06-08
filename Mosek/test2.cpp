#include <bits/stdc++.h>
#include "matrixoperations.h"
using namespace std;

int main() {
    vector<vector< double >> a = {{1,0,-1}};
    vector<vector< double >> b = transpose(a);

    cout <<  "Matrix 1: " << endl;
    for(auto row : a) {
        for(auto col : row) {
            cout << col << " ";
        }
        cout << endl;
    }

    cout << "Matrix 2: " << endl;
    for(auto row : b) {
        for(auto col : row) {
            cout << col << " ";
        }
        cout << endl;
    }

    vector<vector < double >> answer = matrixMult(b,a);
    cout << "Matrix Multiplication result is: " << endl; // (e_ij) * (e_ij^T)

    for(auto row : answer) {
        for(auto col : row) {
            cout << col << " ";
        }
        cout << endl;
    }
}