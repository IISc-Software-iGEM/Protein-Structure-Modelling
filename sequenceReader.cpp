#include <iostream>
#include <fstream>
#include <vector>
#include <vector>
#include <algorithm>
#include <sstream>

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
    // Variable Initialization

    string proteinName = "1qhk";
    string inSeqFile = "1qhk_23_29.seq";
    string proteinPath = proteinName + '/' + inSeqFile;

    cout << "Sequence File Path is: " << proteinPath << endl;

    // Function Call
    pair<string, vector<int>> result = sequenceReader(proteinPath);
    string sequence = result.first;
    vector <int> numberList = result.second;


    cout << sequence << endl;
    for (auto x : numberList) {
        cout << x << " ";
    }
    cout << endl;
}