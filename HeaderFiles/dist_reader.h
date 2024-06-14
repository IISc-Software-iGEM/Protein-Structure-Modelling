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

vector<vector< distanceConstraint* >> constructClassObject(vector<vector<string>> readedUplFiles) {
    vector<vector< distanceConstraint* >> dconstraints;
    for(auto readedFile : readedUplFiles) {

        vector< distanceConstraint* > DCONSTRAINT(readedFile.size());
        for(string x: readedFile) {

            // variable initializations
            int number1, number2;
            string string1, string2, string3, string4;
            float distance, peak;

            // istringstream object to read the line
            istringstream iss(x);

            // reading the line
            iss >> number1 >> string1 >> string2 >> number2 >> string3 >> string4 >> distance;

            // check for peak value
            if(iss) {
                iss >> peak;
            }
            else {
                peak = -1;
            }

            distanceConstraint* distBound = new distanceConstraint(number1, string1, string2, number2, string3, string4, distance, peak);

            DCONSTRAINT.push_back(distBound);
        }

        dconstraints.push_back(DCONSTRAINT);
    }
    // auto X = dconstraints[0];
    // for(auto x : X) {
    //     cout << x->sres << " " << x->stype << " " << x->satom << " " << x->tres << " " << x->ttype << " " << x->tatom << " " << x->distance << " " << x->peak << endl;
    // }
    return dconstraints;
}

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
    //printing the readed lines 
    // for(auto line : lines) {
    //     for(auto l : line) {
    //         cout << l << endl;
    //     }
    // }

    return constructClassObject(lines);
}
