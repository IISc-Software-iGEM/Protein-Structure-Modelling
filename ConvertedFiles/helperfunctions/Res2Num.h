#include <bits/stdc++.h>
using namespace std;

int Res2Num(const char name) {

    vector< char > amino1Letter = {'A','R','N','D','C','E','Q','G','H','I','L','K','M','F','P','S','T','W','Y','V','Z'};

    for(int i = 0; i < amino1Letter.size(); i++) {
        if(name == amino1Letter[i]) {
            return i+1;
        }
    }
    
    return -1;
}