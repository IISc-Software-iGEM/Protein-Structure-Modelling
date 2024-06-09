#include "matrixoperations.h"
#include "boundreader.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

using namespace std;

int main() {
    int dimension = 6;
    int numConstraints = 1;
    
    string equalityBoundFile = "testfile.txt";
    vector < boundReader* > EqualityBounds = bound_reader(equalityBoundFile);

    int n = EqualityBounds.size();
    cout << "Total number of equality bounds are: " << n << endl; 

    Model::t M = new Model("ayushRaina");
    auto _M = finally([&]() {M->dispose();});

    Variable::t XX = M->variable(Domain::inPSDCone(dimension));
    auto traceX = Expr::sum(XX->diag());

    //  Set Objective Function
    M->objective(ObjectiveSense::Minimize, traceX);
    
    for(int i = 0;i < n; i++) {
        boundReader* X = EqualityBounds[i];

        // create E
        vector<vector< double >> E(1, vector<double>(dimension,0));
        E[0][X->x_i-1] = 1;
        E[0][X->x_j-1] = -1;

        // E transpose
        vector<vector <double> > ETranspose = transpose(E);

        // matmul
        vector<vector<double>> answer = matrixMult(ETranspose, E);

        int row = answer.size();
        int col = answer[0].size();

        shared_ptr<ndarray<double, 2>> A = new_array_ptr<double, 2>(shape(row,col));
        for(int i = 0;i < n; i++) { 
            for(int j = 0;j < n; j++) {
                A->operator()(i, j) = answer[i][j];
            }
        }

        auto EZ = Expr::mulDiag(A,XX);
        auto traceEZ = Expr::sum(EZ);

        M->constraint(traceEZ, Domain::lessThan(X->bound));
    }

    M->solve();
    cout << "Solution : " << endl;
    cout << " X = " << *(XX->level()) << endl;
}