#include <bits/stdc++.h>
#include "boundreader.h"
using namespace std;

int main() {
    string upperboundfilename = "upBounds.txt";
    vector < boundReader* > UpperBounds = bound_reader(upperboundfilename);

    // check print
    cout << "Printing the Upper Bounds " << endl;
    for(const auto& X : UpperBounds) {
        cout << X->x_i << " " << X->x_j << " " << X->bound << endl;
    }

    string equalityboundsfilename = "eqBounds.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityboundsfilename);

    // check print
    cout << "Printing the Equality Bounds " << endl;
    for(const auto& X : EqualityBounds) {
        cout << X->x_i << " " << X->x_j << " " << X->bound << endl;
    }

    
}