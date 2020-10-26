#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include "../include/linear_algebra.hpp"

using std::vector, std::function;

// Define this in separate source and load as .so?
double f2(double x, double y) {
//    return sin(y + 1) - x - 1.2;
    return cos(x) + y - 1.5;
}

double f1(double x, double y) {
//    return 2 * y + cos(x) - 2;
    return 2 * x - sin(y - 0.5) - 1;
}

double df1_x(double x, double y) {
//    return 0;
    return sin(x);
}

double df1_y(double x, double y) {
//    return cos(y + 1);
    return 0;
}

double df2_x(double x, double y) {
//    return sin(x) / 2;
    return 0;
}

double df2_y(double x, double y) {
//    return 0;
    return cos(0.5 - y) / 2;
}

double f2N(double x) {
//    return sin(y + 1) - 1.2;
    return 1.5 - cos(x);
}

double f1N(double y) {
//    return 1 - cos(x) / 2;
    return (1 - sin(0.5 - y)) / 2;
}

int main() {
    // Define start values
    double x0 = 0.5;
    double y0 = 0.6;
    double alpha = 1;
    double lambda = 0.5;

    // Read the equation system to solve
    vector<function<double(double, double)>> functions = {f1, f2};
    vector<function<double(double)>> nFunctions(2);
    vector<vector<function<double(double, double)>>> derivatives(2,
                                                                 vector<function<double(double, double)>>(2));
    // Holy crap
    derivatives[0][0] = df1_x;
    derivatives[0][1] = df1_y;
    derivatives[1][0] = df2_x;
    derivatives[1][1] = df2_y;

    nFunctions[0] = f1N;
    nFunctions[1] = f2N;

    // Solve using simple iteration method
    auto res = LAlgebra::NLSimpleIterSolve(x0, y0, functions, nFunctions, derivatives);
    std::cout << res << "\n";

    // Solve using Newton's method

    // Solve using gradient descent

    return 0;
}
