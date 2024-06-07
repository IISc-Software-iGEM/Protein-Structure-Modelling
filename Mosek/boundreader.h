#include <bits/stdc++.h>
#include <algorithm>
using namespace std;

class boundReader {
public:
    double x_i;
    double x_j;
    double bound;
    int flag;

    boundReader(double x, double y, double bound, int zeroOne) {
        this->x_i = x;
        this->x_j = y;
        this->bound = bound;
        this->flag = zeroOne;
    }
};

vector < boundReader* > bound_reader(const string& boundfilename) {
    ifstream file(boundfilename);

    // vector to return
    vector < boundReader* > bounds;

    // check if the file is correctly opened;
    if(!file.is_open()) {
        cerr << "Failed to open the file" << endl;
        return bounds;
    }

    // reading line by line
    string Line;
    while(getline(file, Line)) {
        istringstream iss(Line);

        // read the required values
        double a,b,c;
        int d;

        if(!(iss >> a >> b >> c >> d)) {
            d = -1;
        }

        boundReader* X = new boundReader(a,b,c,d);
        bounds.push_back(X);
    }

    file.close();

    return bounds;
}