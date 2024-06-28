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
    vector< Point* > coordinates;
    vector<int> indexes;

    // function to add new coordinate
    void addCoordinate(Point* P) {
        this->coordinates.push_back(P);
    }

    // function to add index
    void addIndex(int index) {
        this->indexes.push_back(index);
    }
};



int main() {
    vector<  CoordinateAndIndex *  > PatchContent;
    int numPatches = 8;

    vector<string> coordinatePaths;
    for(int i = 0;i < numPatches; i++) {
        string coordPath = "../Groups/group" + to_string(i+1) + ".txt";
        coordinatePaths.push_back(coordPath);
    }

    vector<string> indexPaths;
    for(int i = 0;i < numPatches; i++) {
        string indexPath = "../Indexes/index" + to_string(i+1) + ".txt";
        indexPaths.push_back(indexPath);
    }

    // print coordinate paths
    cout << "Coordinate Paths" << endl;
    for(auto path : coordinatePaths) {
        cout << path << endl;
    }

    // print index paths
    cout << "Index Paths" << endl;
    for(auto path : indexPaths) {
        cout << path << endl;
    }


    // reading the coordinates and indexes
    for(int i = 0;i < numPatches; i++) {
        CoordinateAndIndex* CAI;

        ifstream coordFile(coordinatePaths[i]);
        ifstream indexFile(indexPaths[i]);

        if(!coordFile.is_open() or !indexFile.is_open()) {
            cout << "Error in opening the file" << endl;
            exit(0);
        }
        else {
            cout << "Files are opened successfully" << endl;
            CAI = new CoordinateAndIndex();

            string line;
            while(getline(coordFile, line)) {
                istringstream ss(line);
                double x, y, z;
                ss >> x >> y >> z;
                Point* P = new Point(x, y, z);
                CAI->addCoordinate(P);
            }

            while(getline(indexFile, line)) {
                istringstream ss(line);
                int index;
                ss >> index;
                CAI->addIndex(index);
            }

            PatchContent.push_back(CAI);
        }

    }

    // printing the readed coordinates and indexes
    for(int i = 0;i < numPatches; i++) {
        cout << i+1 << "th Patch" << endl;
        cout << "Number of Coordinates: " << PatchContent[i]->coordinates.size() << endl;
        cout << "Number of Indexes: " << PatchContent[i]->indexes.size() << endl;
    }
}