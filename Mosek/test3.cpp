#include <bits/stdc++.h>
#include "boundreader.h"
#include "matrixoperations.h"

using namespace std;

int main() {
    string equalityboundsfile = "eqbounds.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityboundsfile);

    int n = EqualityBounds.size();

    cout << "Total number of equality bounds are: " << n << endl;

    for(int i = 0;i < n; i++) {
        boundReader* X = EqualityBounds[i];

        // Create E
        vector<vector<double>> E(1, vector<double>(n,0));
        E[0][X->x_i-1] = 1;
        E[0][X->x_j-1] = -1;

        // E transpose
        vector< vector <double> > ETranspose = transpose(E);
        
        // Matmul
        vector<vector<double>> answer = matrixMult(ETranspose,E);
        

    }
}
