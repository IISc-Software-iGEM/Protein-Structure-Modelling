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

void displayCoordinates( CoordinateAndIndex* CAI) {
    for(int i = 0;i < CAI->coordinates.size(); i++) {
        cout << CAI->coordinates[i]->x << " " << CAI->coordinates[i]->y << " " << CAI->coordinates[i]->z << endl;
    }
}

void displayIndexes( CoordinateAndIndex* CAI) {
    for(int i = 0;i < CAI->indexes.size(); i++) {
        cout << CAI->indexes[i] << " ";
    }
}

void printPaths(vector<string> paths) {
    for(auto path : paths) {
        cout << path << endl;
    }
}

vector< CoordinateAndIndex* > constructTheStruct(const vector<string> coordinatePaths, const vector<string> indexPaths) {
    vector< CoordinateAndIndex* > PatchContent;
    int numPatches = coordinatePaths.size();

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

    // uncomment to print the content of kth Patch
    // int k = 3;
    // displayCoordinates(PatchContent[k-1]);
    // displayIndexes(PatchContent[k-1]);
    return PatchContent;
}

vector<int> mapToTmpIndex(const vector <int> globalIndex, const vector<int> atom_mapVector) {
    int n_atom = globalIndex.size();
    vector<int> tmpIndex(n_atom,0);

    for(int i = 0;i < n_atom; i++) {
        int index = globalIndex[i];
        
        auto forwardIterator = find(atom_mapVector.begin(), atom_mapVector.end(), index);
        
        if(forwardIterator != atom_mapVector.end()) {
            tmpIndex[i] = distance(atom_mapVector.begin(), forwardIterator) + 1;
        }
        else {
            cout << "Error in mapping the indexes" << endl;
        }
    }

    return tmpIndex;
}

pair <vector < CoordinateAndIndex* >, int> doRegForGroup(const vector<int> groupNums, const vector < CoordinateAndIndex* > CAI) {
    
    if(groupNums.size() != CAI.size()) {
        cout << "Number of groups and patches are not equal" << endl;
        exit(0);
    }
    else {
        cout << "Number of groups and patches are equal" << endl;
    }

    set<int> atom_map;
    for(int i = 0;i < groupNums.size(); i++) {
        vector<int> Indexes = CAI[i]->indexes;   
        for(int j = 0;j < Indexes.size(); j++) {
            atom_map.insert(Indexes[j]);
        }
    }

    // convert set into a vector
    vector<int> atom_mapVector(atom_map.begin(), atom_map.end());
    atom_map.clear();

    // Uncomment to print the contents of atom_mapVector
    // for(auto x : atom_mapVector) {
    //     cout << x << " ";
    // }
    // cout << endl;
    cout << "Unique Elements: " << atom_mapVector.size() << endl;
    
    // mapping the global indexes to temporary indexes
    // auto mappedIndex = mapToTmpIndex(CAI[5]->indexes, atom_mapVector);
    // cout << mappedIndex.size() << endl;
    

    // Map to TMP index for all the patches
    for(int i = 0;i < CAI.size(); i++) {
        CAI[i]->indexes = mapToTmpIndex(CAI[i]->indexes, atom_mapVector);
    }

    return {CAI, atom_mapVector.size()};
}

void formBL(vector< CoordinateAndIndex* > CAI, int dimension, int uniqueIndex) {
    
    int m = CAI.size();
    int Md = m*dimension;

    vector<vector<double>> bTemp(uniqueIndex+m, vector<double>(dimension, 0));
    cout << "Rows: " << uniqueIndex+m << " Columns: " << dimension << endl;

    vector< Point* > Coordinates = CAI[1]->coordinates;
    vector<int> Indexes = CAI[1]->indexes;

    for(int i = 0;i < Indexes.size(); i++) {
        bTemp[Indexes[i]-1][0] = Coordinates[i]->x;
        bTemp[Indexes[i]-1][1] = Coordinates[i]->y;
        bTemp[Indexes[i]-1][2] = Coordinates[i]->z;
    }

    for(auto row : bTemp) {
        for(auto x : row) {
            cout << x << " ";
        }
        cout << endl;
    }

    return;
}

void doGretSDP(vector < CoordinateAndIndex* > CAI, int dimension, int uniqueIndexes) {
    formBL(CAI, dimension, uniqueIndexes);
    return;
}

int main() {
    int dimension = 3;

    vector<string> coordinatePaths;
    for(int i = 0;i < 8; i++) {
        string coordPath = "../Groups/group" + to_string(i+1) + ".txt";
        coordinatePaths.push_back(coordPath);
    }

    vector<string> indexPaths;
    for(int i = 0;i < 8; i++) {
        string indexPath = "../Indexes/index" + to_string(i+1) + ".txt";
        indexPaths.push_back(indexPath);
    }
    
    // construct the struct
    vector< CoordinateAndIndex* > PatchContent = constructTheStruct(coordinatePaths, indexPaths);
    vector<int> groupNums = {1,2,3,4,6,7,42,63};

    // map to temp index
    auto result = doRegForGroup(groupNums, PatchContent);
    PatchContent = result.first;
    int uniqueIndexes = result.second;

    // call the Global Registration module
    doGretSDP(PatchContent, dimension, uniqueIndexes);
}