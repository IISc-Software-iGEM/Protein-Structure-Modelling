#include <bits/stdc++.h>
using namespace std;

class Point {
public:
    double x, y, z;

    // constructor
    Point(double X, double Y, double Z) {
        this->x = X;
        this->y = Y;
        this->z = Z;
    }
};

class CoordinateAndIndex {
public:
    int sizeOfPatch;
    vector<vector<double>> coordinates;
    vector<int> indexes;

    // constructor
    CoordinateAndIndex(int PatchSize) {
        this->sizeOfPatch = PatchSize;
        this->coordinates.resize(PatchSize, vector<double>(3, 0));
        this->indexes.resize(PatchSize, 0);
    }

    // function to add new coordinate
    void addCoordinate(Point* P, int indexNumber) {
        this->coordinates[indexNumber][0] = P->x;
        this->coordinates[indexNumber][1] = P->y;
        this->coordinates[indexNumber][2] = P->z;
    }

    // function to add index
    void addIndex(int index, int indexNumber) {
        this->indexes[indexNumber] = index;
    }
};



int main() {
    
}