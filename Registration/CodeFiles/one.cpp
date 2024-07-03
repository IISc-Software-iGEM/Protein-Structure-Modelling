#include <bits/stdc++.h>
#include "../../Mosek/matrixoperations.h"
#include "../../../eigen-3.4.0/Eigen/Dense"

using namespace std;
using namespace Eigen;


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

vector<vector<double>> formB_matrix(vector< CoordinateAndIndex* > CAI, int dimension, int uniqueIndex) {
    
    int m = CAI.size();
    int Md = m*dimension;
    
    // initialize B matrix with 0 
    vector<vector<double>> B(uniqueIndex+m, vector<double>(Md, 0));

    for(int i = 0;i < m;i++) {
        vector < Point* > Coordinates = CAI[i]->coordinates;
        vector < int > Indexes = CAI[i]->indexes;

        for(int j = 0;j < Indexes.size(); j++) {
            B[Indexes[j]-1][i*dimension] = Coordinates[j]->x;
            B[Indexes[j]-1][i*dimension+1] = Coordinates[j]->y;
            B[Indexes[j]-1][i*dimension+2] = Coordinates[j]->z;
        }
    }

    return transpose(B);
}

pair < vector<vector<double>>, vector<vector<double>> > form_L_and_Adj_matrix(vector< CoordinateAndIndex* > CAI, int dimension, int uniqueIndex) {
    
    int m = CAI.size();

    // initialize L and ADJ matrix with 0
    vector<vector<double>> L(uniqueIndex+m, vector<double>(uniqueIndex+m,0));
    vector<vector<double>> Adj(uniqueIndex+m, vector<double>(uniqueIndex+m, 0));

    for(int i = 0;i < m; i++) {
        vector< int > Indexes = CAI[i]->indexes;

        for(int j = 0;j < Indexes.size(); j++) {
            int currIndex = Indexes[j];

            // filling L matrix
            L[currIndex-1][currIndex-1] += 1;
            L[uniqueIndex+i][uniqueIndex+i] += 1;
            L[currIndex-1][uniqueIndex+i] += -1;
            L[uniqueIndex+i][currIndex-1] += -1;

            // filling Adj matrix
            Adj[currIndex-1][uniqueIndex+i] = 1;
            Adj[uniqueIndex+i][currIndex-1] = 1;
        }

    }
    return {L, Adj};
}

MatrixXd convertToEigenMatrix(vector<vector<double>> matrix) {

    int rows = matrix.size();
    int cols = matrix[0].size();

    MatrixXd M(rows, cols);
    for(int i = 0;i < rows; i++) {
        for(int j = 0;j < cols; j++) {
            M(i,j) = matrix[i][j];
        }
    }

    return M;
}

MatrixXd computePseudoinverse(const MatrixXd &matrix) {
    // Compute the SVD of the input matrix
    JacobiSVD<MatrixXd> svd(matrix, ComputeThinU | ComputeThinV);
    
    // Get the singular values and the U, V matrices
    const VectorXd &singularValues = svd.singularValues();
    const MatrixXd &U = svd.matrixU();
    const MatrixXd &V = svd.matrixV();

    // Create a diagonal matrix for the inverse of singular values
    VectorXd singularValuesInv(singularValues.size());
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {  // Use a threshold to avoid division by zero
            singularValuesInv(i) = 1.0 / singularValues(i);
        } else {
            singularValuesInv(i) = 0;
        }
    }

    // Compute the pseudoinverse using the SVD components
    MatrixXd pseudoinverse = V * singularValuesInv.asDiagonal() * U.transpose();
    return pseudoinverse;
}

bool isSymmetric(const MatrixXd &matrix) {
    return matrix.isApprox(matrix.transpose());
}

void doGretSDP(vector < CoordinateAndIndex* > CAI, int dimension, int uniqueIndexes) {
    vector<vector<double>> B = formB_matrix(CAI, dimension, uniqueIndexes);
    cout << "B Matrix is formed" << endl;
    cout << "Rows: " << B.size() << " and Columns: " << B[0].size() << endl;

    auto pairBL = form_L_and_Adj_matrix(CAI, dimension, uniqueIndexes);
    vector<vector<double>> L = pairBL.first;
    vector<vector<double>> Adj = pairBL.second;

    cout << "L and Adj Matrix is formed" << endl;
    
    // converting B, L and Adj matrix to Eigen Matrix
    MatrixXd B_eigen = convertToEigenMatrix(B);
    MatrixXd L_eigen = convertToEigenMatrix(L);
    MatrixXd Adj_eigen = convertToEigenMatrix(Adj);

    cout << "Matrices converted to Eigen Matrices successfully " << endl;
    cout << "Dimensions of B: " << B_eigen.rows() << " " << B_eigen.cols() << endl;
    cout << "Dimensions of L: " << L_eigen.rows() << " " << L_eigen.cols() << endl;
    cout << "Dimensions of Adj: " << Adj_eigen.rows() << " " << Adj_eigen.cols() << endl;


    MatrixXd B_Ldiag_BT = B_eigen * computePseudoinverse(L_eigen) * B_eigen.transpose();
    
    // print first ith row of B_Ldiag_BT
    cout << "B_Ldiag_BT is computed successfully" << endl;
    cout << "Dimensions of B_Ldiag_BT: " << B_Ldiag_BT.rows() << " " << B_Ldiag_BT.cols() << endl;

    // check is matrix formed is symmetric or not
    if(!isSymmetric(B_Ldiag_BT)) {
        cout << "Matrix is not symmetric" << endl;
        cout << "Making the matrix symmetric" << endl;
        B_Ldiag_BT = 0.5 * (B_Ldiag_BT + B_Ldiag_BT.transpose());
    }

    // for(int i = 0;i < 24; i++) {
    //     cout << "Col " << i+1 << " --> " << B_Ldiag_BT(9,i) << endl;
    // }
    cout << endl;

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