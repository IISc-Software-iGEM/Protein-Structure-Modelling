#include <bits/stdc++.h>
#include <algorithm>
#include "ConvertedFiles/inputreaders/Hbond.h"
#include "ConvertedFiles/inputreaders/dist_reader.h"
using namespace std;

pair<vector<int>, vector<int>> sequenceReader(const string &filePath) {
    // ifstream file(filePath);
    // if(!file.is_open()) {
    //     cerr << "Unable to open file: " << filePath << endl;
    //     return {"", {}};
    // }

    // string sequence;
    // vector<int> numberList;
    // string line;

    // while(getline(file,line)) {
    //     istringstream iss(line);
    //     string word;
    //     int number;

    //     while(iss >> word >> number) {
    //         sequence += word;
    //         sequence += " ";
    //         numberList.push_back(number);
    //     }
    // }

    // file.close();
    // return {sequence,numberList};

    vector<int> sequence = {18,3,6,5,12,7,3};
    vector<int> numberList = {23,24,25,26,27,28,29};


    return {sequence,numberList};
}

int main() {

    // Drop side-chain hydrogen atoms or not
    
    //int hydrogenOmission = 1;
    int hydrogenOmission = 0;
    
    // Number of BFGS iterations
    double gdTolerance = pow(10,-9);
    int cgIterations = 500;
    vector<int> f = {10,10,10,10,10}; /* f_hb, f_tau, f_tal ,f_vdw, f_tas */

    // Output
    cout << "=========================================" << endl;
    cout << "SDP-Based Protein Structure Determination" << endl;

    // Start the timer
    auto timeStart = chrono::high_resolution_clock::now();

    // Variable Initializations
    int inMaxRes;
    int inMinRes;

    string proteinName = "1qhk";
    string inSeqFile = "1qhk_23_29.seq";
    vector<string> inUplFile = {"1qhk_concat_dist_23_29.upl"};
    string inHbondFile = "1qhk_concat_hBond_23_29.upl";

    // Printing to console
    cout << "**********************************************************" << endl;
    cout << "Protein: " << proteinName << endl;
    cout << "Reading input files ------- " << endl;

    string ProteinPath = proteinName + '/' + inSeqFile;
    pair<vector<int>, vector<int>> result = sequenceReader(ProteinPath);
    vector<int> sequence = result.first;
    vector <int> numberList = result.second;

    // Printing the sequence
    cout << "sequence = " << endl;
    for(auto element : sequence) {
        cout << element << " ";
    }
    cout << endl;

    // Printing the number list
    cout << "numberList = " << endl;
    for(auto element : numberList) {
        cout << element << " ";
    }
    cout << endl;

    inMaxRes = *max_element(numberList.begin(), numberList.end());
    inMinRes = *min_element(numberList.begin(), numberList.end());

    cout << inMaxRes << " " << inMinRes << endl;

    // Number of Upl Files
    size_t numUplFiles = inUplFile.size(); 
    cout << "Number of UPL files are >>>> " << inUplFile.size() << endl;

    // Vector for UPL File paths
    vector<string> uplFilePaths(numUplFiles);

    // Adding paths
    for(size_t i = 0;i < numUplFiles; i++) {
        uplFilePaths[i] = proteinName + '/' + inUplFile[i];
    }

    // Printing the upl file paths
    cout << "Paths for UPL file >>>>>> " << endl;
    for(auto &file : uplFilePaths) {
        cout << file << endl;
    }

    if(!inHbondFile.empty()) {
        string hbondFile = proteinName + '/' + inHbondFile;
        string hbondWriteFile = proteinName + "/" + proteinName + "_hbo.upl";

        cout << "Hydrogen Bond given and write files >>>>>>> " << endl;
        cout << hbondFile << endl << hbondWriteFile << endl;

        hbond_reader_me(hbondFile, hbondWriteFile);
        numUplFiles += 1;
        uplFilePaths.push_back(hbondWriteFile);
    }

    cout << "***************" << endl;
    cout << "Updated UPL file paths >>>>>>> " << endl;
    for(auto &file : uplFilePaths) {
        cout << file << endl;
    }

    // Reading the UPL files
    vector<vector< distanceConstraint* >> raw_up = readUplFiles(uplFilePaths);

    // Printing the raw_up
    cout << "***************" << endl;
    cout << "Printing the raw_up >>>>>>> " << endl;
    cout << raw_up.size() << endl; // 74 struct and 2 struct.

    // check for correct working
    // for(auto up : raw_up) {
    //     for(auto up2 : up) {
    //         cout << up2 -> sres << " " << up2 -> stype << " " << up2 -> satom << " " << up2 -> tres << " " << up2 -> ttype << " " << up2 -> tatom << " " << up2 -> distance << " " << up2 -> peak << endl;
    //     }
    // }




    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(end - timeStart);
    cout << "Time taken to read input files: " << duration.count() << " seconds" << endl;
    return 0;
}
