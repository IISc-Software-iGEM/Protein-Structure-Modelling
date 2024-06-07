#include <iostream>
#include <random>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;


int main() { 
    int dimension = 133;
    int numConstraints = 249;
    
    // we have 76 upper bounds and 173 equality bounds + 1 lower bound

    Model::t M = new Model("ayushRaina");
    auto _M = finally([&]() {M->dispose();});

    Variable::t X = M->variable(Domain::inPSDCone(dimension));
    auto traceX = Expr::sum(X->diag());

    //  Set Objective Function
    M->objective(ObjectiveSense::Minimize, traceX);

    

}