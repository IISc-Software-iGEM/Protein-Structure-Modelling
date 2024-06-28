#include <iostream>
#include <vector>
#include <string>

char Three2One(const std::string& input) {
    std::vector<std::string> amino_1L = {"A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "C"};
    std::vector<std::string> amino_3L = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "CYS"};
    char type = ' ';
    for (int i = 0; i < amino_3L.size(); i++) {
        if (amino_3L[i] == input) {
            type = amino_1L[i][0];
            break;
        }
    }
    return type;
}
