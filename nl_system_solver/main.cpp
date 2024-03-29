#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include "../include/linear_algebra.hpp"

using std::vector, std::function;

// Define this in separate source and load as .so?
double f1(double x, double y) {
//    return sin(y - 1) + x - 1.3;
    return cos(x) + y - 1.5;
}

double f2(double x, double y) {
//    return y - sin(x + 1) - 0.8;
    return 2 * x - sin(y - 0.5) - 1;
}

double df1_x(double x, double y) {
    return sin(x);
//    return 0;
}

double df1_y(double x, double y) {
    return 0;
//    return -cos(y - 1);
}

double df2_x(double x, double y) {
    return 0;
//    return cos(x + 1);
}

double df2_y(double x, double y) {
    return cos(0.5 - y) / 2;
//    return 0;
}

double ddf1_x(double x, double y) {
    return -sin(x);
//    return 1;
}

double ddf1_y(double x, double y) {
    return 1;
//    return cos(y - 1);
}

double ddf2_x(double x, double y) {
    return 2;
//    return -cos(x + 1);
}

double ddf2_y(double x, double y) {
    return -cos(0.5 - y);
//    return 1;
}

double ddf1(double x, double y) {
    return 4 * (sin(0.5 - y) + 2 * x - 1) - 2 * sin(x) * (cos(x) + y - 1.5);
//    return 2 * (-cos(x + 1) * (-sin(x + 1) + y - 0.8) + x - sin(1 - y) - 1.3);
}

double ddf2(double x, double y) {
    return 2 * ((-cos(0.5 - y)) * (sin(0.5 - y) + 2 * x - 1) + cos(x) + y - 1.5);
//    return 2 * (cos(1 - y) * (x - sin(1 - y) - 1.3) - sin(x + 1) + y - 0.8);
}

double f1N(double y) {
    return (1 - sin(0.5 - y)) / 2;
//    return 1.3 - sin(y - 1);
}

double f2N(double x) {
//    return 0.8 + sin(x + 1);
    return 1.5 - cos(x);
}


int main() {
    // Define start values
    double x0 = 0.6;
    double y0 = 0.7;
    double alpha = 1;
    double lambda = 0.5;

    // Read the equation system to solve
    vector<function<double(double, double)>> functions = {f1, f2};
    vector<function<double(double)>> nFunctions(2);
    vector<vector<function<double(double, double)>>> derivatives(2,
                                                                 vector<function<double(double, double)>>(2));
    vector<vector<function<double(double, double)>>> derivatives1(2,
                                                                  vector<function<double(double, double)>>(2));
    vector<function<double(double, double)>> derivatives2(2);

    // Holy crap
    derivatives[0][0] = df1_x;
    derivatives[0][1] = df1_y;
    derivatives[1][0] = df2_x;
    derivatives[1][1] = df2_y;

    derivatives1[0][0] = ddf1_x;
    derivatives1[0][1] = ddf1_y;
    derivatives1[1][0] = ddf2_x;
    derivatives1[1][1] = ddf2_y;

    derivatives2[0] = ddf1;
    derivatives2[1] = ddf2;

    nFunctions[0] = f1N;
    nFunctions[1] = f2N;

    std::cout << "Start at " << std::setprecision(4) << x0 << " " << y0 << "\n";

    // Solve using simple iteration method
    auto res = LAlgebra::NLSimpleIterSolve(x0, y0, functions, nFunctions, derivatives);

    // Solve using Newton's method
    std::cout << "\nNewton's method:\n";
    res = LAlgebra::NLNewtonSolve(x0, y0, functions, derivatives, derivatives1);

    // Solve using gradient descent
    std::cout << "\nGradient descent:\n";
    res = LAlgebra::NLGradientDescentSolve(x0, y0, alpha, lambda, functions, derivatives2, derivatives);

    return 0;
}
