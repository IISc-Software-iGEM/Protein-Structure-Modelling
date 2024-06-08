#include <bits/stdc++.h>
#include "boundreader.h"
#include "matrixoperations.h"

using namespace std;

int main() {
    string equalityboundsfile = "testfile.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityboundsfile);

    int n = EqualityBounds.size();

    cout << "Total number of equality bounds are: " << n << endl;

    vector<int> rows, cols;
    vector<double> values;

    for(int i = 0;i < n; i++) {
        boundReader* X = EqualityBounds[i];

        // Create E
        vector<vector<double>> E(1, vector<double>(n,0));
        E[0][X->x_i-1] = 1;
        E[0][X->x_j-1] = -1;

        // for(auto row : E) {
        //     for(auto col : row) {
        //         cout << col << " ";
        //     }
        //     cout << endl;
        // }
        
        // E transpose
        vector< vector <double> > ETranspose = transpose(E);
        // for(auto row : ETranspose) {
        //     for(auto col : row) {
        //         cout << col << " ";
        //     }
        //     cout << endl;
        // }

        // Matmul
        vector<vector<double>> answer = matrixMult(ETranspose,E);
        for(auto row : answer) {
            for(auto col : row) {
                cout << col << " ";
            }
            cout << endl;
        }

        // we will create sparse matrices
        cout << "row vector for mosek: " << endl;
        rows.push_back(X->x_i);
        rows.push_back(X->x_j);
        rows.push_back(X->x_i);
        rows.push_back(X->x_j);
        for(auto x : rows) {
            cout << x << " ";
        }
        cout << endl;

        cout << "column vector for mosek: " << endl;
        cols.push_back(X->x_i);
        cols.push_back(X->x_j);
        cols.push_back(X->x_j);
        cols.push_back(X->x_i);
        for(auto x : cols) {
            cout << x << " ";
        }
        cout << endl;

        cout << "values vector for mosek: " << endl;
        values.push_back(1);
        values.push_back(1);
        values.push_back(-1);
        values.push_back(-1);
        for(auto x : values) {
            cout << x << " ";
        }
        cout << endl;
        
        cout << " ============================== " << endl;

        rows.clear();
        cols.clear();
        values.clear();
    }
}
