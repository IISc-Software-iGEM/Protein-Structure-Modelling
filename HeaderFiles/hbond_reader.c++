#include <bits/stdc++.h>
using namespace std;

class HBondData {
    int value1;
    string str1;
    string str2;

    int value2;
    string str3;
    string str4;

    int value3;
};

void hbond_reader_me(const string& fileA) {
    throw invalid_argument("Function requires exactly 2 arguments");
}

void hbond_reader_me(const string& fileA, const string& fileB, const string& fileC) {
    throw invalid_argument("Function requires exactly 2 arguments");
}

void hbond_reader_me() {
    throw invalid_argument("Function requires exactly 2 arguments");
}

void hbond_reader_me(const string& fileName, const string& writeFileName) {

    // this function should work
    cout << "Hydrogen Bond Reader called properly" << endl;

    // open the file and throw exception if file could not be opened.
    ifstream infile(fileName);
    if(!infile) {
        throw runtime_error("Hbond file could not be opened");
    }

    // create vector to store the raw data
    vector < HBondData > rawData;

    // parsing the file
    string line;
    while(getline(infile, line)) {
        istringstream iss(line);
        HBondData data;
        if()
    }
}