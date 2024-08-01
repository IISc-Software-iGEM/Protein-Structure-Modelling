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

void timeDisplay(chrono::high_resolution_clock::time_point startTime, chrono::high_resolution_clock::time_point endTime, string message, int i) {
    cout << "Time Display " << i << ": " << endl;
    cout << "Time taken for " << message << ": " << chrono::duration_cast<chrono::microseconds>(endTime-startTime).count() << " microseconds" << endl;
}

void totalDisplay(chrono::high_resolution_clock::time_point startTime, chrono::high_resolution_clock::time_point endTime) {
    cout << "Total Time: " << chrono::duration_cast<chrono::microseconds>(endTime-startTime).count() << " microseconds" << endl;
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
            // cout << "Error in opening the file" << endl;
            exit(0);
        }
        else {
            // cout << "Files are opened successfully" << endl;
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
    // for(int i = 0;i < numPatches; i++) {
    //     cout << i+1 << "th Patch" << endl;
    //     cout << "Number of Coordinates: " << PatchContent[i]->coordinates.size() << endl;
    //     cout << "Number of Indexes: " << PatchContent[i]->indexes.size() << endl;
    // }

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
        if (singularValues(i) > 1e-6) {  // Use a threshold to avoid division by zero
            singularValuesInv(i) = 1.0 / singularValues(i);
        } else {
            singularValuesInv(i) = 0;
        }
    }

    // Compute the pseudoinverse using the SVD components
    MatrixXd pseudoinverse = V * singularValuesInv.asDiagonal() * U.transpose();

    auto endTime2 = chrono::high_resolution_clock::now();
    return pseudoinverse;
}

bool isSymmetric(const MatrixXd &matrix) {
    return matrix.isApprox(matrix.transpose());
}

void saveMatrix(const MatrixXd& matrix, const std::string& filename) {
    std::ofstream file(filename);
    if (file.is_open()) {
        file << matrix.rows() << " " << matrix.cols() << "\n";
        file << matrix << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
    }
}

void doGretSDP(vector < CoordinateAndIndex* > CAI, int dimension, int uniqueIndexes) {

    auto startTime = chrono::high_resolution_clock::now();
    vector<vector<double>> B = formB_matrix(CAI, dimension, uniqueIndexes);
    cout << "B Matrix is formed" << endl;
    cout << "Rows: " << B.size() << " and Columns: " << B[0].size() << endl;
    
    
    auto pairBL = form_L_and_Adj_matrix(CAI, dimension, uniqueIndexes);
    vector<vector<double>> L = pairBL.first;
    vector<vector<double>> Adj = pairBL.second;

    auto endTime2 = chrono::high_resolution_clock::now();
    timeDisplay(startTime, endTime2, "forming B, L and Adj matrix", 2);
    // cout << "L and Adj Matrix is formed" << endl;
    
    // // converting B, L and Adj matrix to Eigen Matrix
    MatrixXd B_eigen = convertToEigenMatrix(B);
    MatrixXd L_eigen = convertToEigenMatrix(L);
    MatrixXd Adj_eigen = convertToEigenMatrix(Adj);

    auto endTime3 = chrono::high_resolution_clock::now();
    timeDisplay(endTime2, endTime3, "converting to Eigen Matrix", 3);

    saveMatrix(L_eigen, "L_eigen.txt");
    saveMatrix(B_eigen, "B_eigen.txt");

    auto endTime4 = chrono::high_resolution_clock::now();
    timeDisplay(endTime3, endTime4, "saving the Eigen Matrix", 4);
    
    cout << "Matrices converted to Eigen Matrices successfully " << endl;
    cout << "Dimensions of B: " << B_eigen.rows() << " " << B_eigen.cols() << endl;
    cout << "Dimensions of L: " << L_eigen.rows() << " " << L_eigen.cols() << endl;
    cout << "Dimensions of Adj: " << Adj_eigen.rows() << " " << Adj_eigen.cols() << endl;

    auto endTime5 = chrono::high_resolution_clock::now();
    timeDisplay(endTime4, endTime5, "displaying the dimensions of Eigen Matrix", 5);

    // computing the pseudoinverse of L
    MatrixXd L_pseudo = computePseudoinverse(L_eigen);

    auto endTime6 = chrono::high_resolution_clock::now();
    timeDisplay(endTime5, endTime6, "computing the pseudoinverse of L", 6);

    MatrixXd B_Ldiag_BT = B_eigen * L_pseudo * B_eigen.transpose();
    auto endTime7 = chrono::high_resolution_clock::now();
    timeDisplay(endTime6, endTime7, "computing B_Ldiag_BT", 7);

    // // print first ith row of B_Ldiag_BT
    // cout << "B_Ldiag_BT is computed successfully" << endl;
    // cout << "Dimensions of B_Ldiag_BT: " << B_Ldiag_BT.rows() << " " << B_Ldiag_BT.cols() << endl;

    // check is matrix formed is symmetric or not
    if(!isSymmetric(B_Ldiag_BT)) {
        cout << "Matrix is not symmetric" << endl;
        cout << "Making the matrix symmetric" << endl;
        B_Ldiag_BT = 0.5 * (B_Ldiag_BT + B_Ldiag_BT.transpose());
    }
    auto endTime8 = chrono::high_resolution_clock::now();
    timeDisplay(endTime7, endTime8, "checking the symmetry of B_Ldiag_BT", 8);

    // // for(int i = 0;i < 24; i++) {
    // //     cout << "Col " << i+1 << " --> " << B_Ldiag_BT(9,i) << endl;
    // // }

    std::ofstream file("objectiveMatrix.txt");
    if (file.is_open()) {
        file << B_Ldiag_BT.rows() << " " << B_Ldiag_BT.cols() << "\n";
        file << B_Ldiag_BT << "\n";
        file.close();
    } else {
        std::cerr << "Unable to open file for writing: " << "objectiveMatrix.txt" << std::endl;
    }
    cout << endl;

    auto endTime9 = chrono::high_resolution_clock::now();
    timeDisplay(endTime8, endTime9, "saving the objective matrix", 8);

    return;
}

int main() {
    int dimension = 3;


    auto startTime = chrono::high_resolution_clock::now();

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
    
    
    auto endTime1 = chrono::high_resolution_clock::now();
    timeDisplay(startTime, endTime1, "writing input file paths", 1);
    totalDisplay(startTime, endTime1);

    // construct the struct
    vector< CoordinateAndIndex* > PatchContent = constructTheStruct(coordinatePaths, indexPaths);
    
     
    auto endTime2 = chrono::high_resolution_clock::now();
    timeDisplay(endTime1, endTime2, "reading the files and constructing the struct", 2);
    totalDisplay(startTime, endTime2);

    vector<int> groupNums = {1,2,3,4,6,7,42,63};

    // map to temp index
    auto result = doRegForGroup(groupNums, PatchContent);
    PatchContent = result.first;
    int uniqueIndexes = result.second;

    auto endTime3 = chrono::high_resolution_clock::now();
    timeDisplay(endTime2, endTime3, "mapping to temporary index", 3);
    totalDisplay(startTime, endTime3);

    // call the Global Registration module
    doGretSDP(PatchContent, dimension, uniqueIndexes);

    auto endTime4 = chrono::high_resolution_clock::now();
    timeDisplay(endTime3, endTime4, "Global Registration", 4);
    totalDisplay(startTime, endTime4);
}