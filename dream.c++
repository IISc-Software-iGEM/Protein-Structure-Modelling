#include <bits/stdc++.h>
#include <algorithm>
using namespace std;

pair<string, vector<int>> sequenceReader(const string &filePath) {
    ifstream file(filePath);
    if(!file.is_open()) {
        cerr << "Unable to open file: " << filePath << endl;
        return {"", {}};
    }

    string sequence;
    vector<int> numberList;
    string line;

    while(getline(file,line)) {
        istringstream iss(line);
        string word;
        int number;

        while(iss >> word >> number) {
            sequence += word;
            sequence += " ";
            numberList.push_back(number);
        }
    }

    file.close();
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
    pair<string, vector<int>> result = sequenceReader(ProteinPath);
    string sequence = result.first;
    vector <int> numberList = result.second;

    inMaxRes = *max_element(numberList.begin(), numberList.end());
    inMinRes = *min_element(numberList.begin(), numberList.end());

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
    }
    return 0;
}
