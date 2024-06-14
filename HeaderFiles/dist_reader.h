#include <bits/stdc++.h>
using namespace std;

class distanceConstraint {
public:

    int sres;
    string stype;
    string satom;

    int tres;
    string ttype;
    string tatom;

    float distance;
    float peak;

    // Constructor
    distanceConstraint(int num1, string type1, string atom1, int num2, string type2, string atom2, float DISTANCE, float PEAK=0) {
        this->sres = num1;
        this->stype = type1;
        this->satom = atom1;

        this->tres = num2;
        this->ttype = type2;
        this->tatom = atom2;

        this->distance = DISTANCE;
        this->peak = PEAK;
    }
};

vector<vector< distanceConstraint* >> readUplFiles(const vector<string> uplFilePaths) {
    
    // check for upl files are received or not
    if(uplFilePaths.size() == 0) {
        cout << "No UPL files are received" << endl;
        exit(0);
    }
    else {
        
        cout << "***************UPL files are received***************" << endl;
        for(auto file : uplFilePaths) {
            cout << file << endl;
        }
    }

    // vector to store the readed lines
    vector<vector<string>> lines;

    cout << "***************New UPL file paths***************" << endl;
    for(auto file : uplFilePaths) {
        cout << "Reading >>>>>> " << file << endl;
        vector<string> LINE;

        ifstream filename(file);
        if(!filename.is_open()) {
            cout << "Error in opening the file >>>>>> " << file << endl;
            exit(0);
        }
        else {
    
            cout << "File: " << file << " is opened successfully" << endl;
            string line;

            while(getline(filename, line)) {
                if(!line.empty() and line.back() == '\r') {
                    line.pop_back();
                }
                LINE.push_back(line);
            }
        }

        lines.push_back(LINE);
    }

    cout << "***************Lines are readed successfully***************" << endl;

    vector<vector< distanceConstraint* >> distCons;

    // printing the readed lines 
    for(auto line : lines) {
        vector< distanceConstraint* > one;
        for(auto l : line) {
            istringstream ss(l);
            int num1, num2;
            string type1, type2, atom1, atom2;
            float distance, peak;
            
            int wordCount = 0;
            istringstream temp(l);
            string word;

            while(temp >> word) {
                wordCount++;
            }

            if(wordCount == 7) {
                ss >> num1 >> type1 >> atom1 >> num2 >> type2 >> atom2 >> distance;
                peak = 0;
            }
            else if(wordCount == 9) {
                ss >> num1 >> type1 >> atom1 >> num2 >> type2 >> atom2 >> distance >> peak;
                peak = -1;
            }
            
            distanceConstraint* dc = new distanceConstraint(num1, type1, atom1, num2, type2, atom2, distance, peak);
            one.push_back(dc);
        }
        distCons.push_back(one);
    }

    // printing the distance constraints
    return distCons;
}
