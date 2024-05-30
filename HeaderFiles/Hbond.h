#include <bits/stdc++.h>
using namespace std;

class HBondData {
public:
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
        
        if(!(iss >> data.value1 >> data.str1 >> data.str2 >> data.value2 >> data.str3 >> data.str4 >> data.value3)) {

        }

        rawData.push_back(data);
    }

    // closing the input file
    infile.close();
    

    // opening the output files
    ofstream outfile(writeFileName);
    if(!outfile) {
        throw runtime_error("Output file could not be opened");
    }

    // process and write the data to output file
    for (const HBondData& data : rawData) {
        outfile << data.value1 << " " << data.str1 << " " << data.str2 << '\t' << data.value2 << " " << data.str3 << " " << data.str4 << data.value3 << "\t\t#peak\t-1\n";
    }

    //closing the output file
    outfile.close();

}