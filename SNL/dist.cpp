#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

struct Point {
    double x, y, z;
};

double calculateDistance(const Point& p1, const Point& p2) {
    return pow(p2.x - p1.x, 2) + 
                pow(p2.y - p1.y, 2) + 
                pow(p2.z - p1.z, 2);
}

int main() {
    vector<Point> points;
    ifstream inputFile("points.txt");
    if (!inputFile) {
        cerr << "Unable to open file" << endl;
        return 1;
    }

    Point p;
    while (inputFile >> p.x >> p.y >> p.z) {
        points.push_back(p);
    }
    inputFile.close();

    int n = points.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double distance = calculateDistance(points[i], points[j]);
            cout << "Distance between point " << i+1 << " and point " << j+1 << " is: " << distance << endl;
        }
    }

    return 0;
}
