#include <cmath>
#include <iostream>
#include <iomanip>
#include "../include/computational_math.hpp"

using std::pow, std::cout;

double func1(double x) {
    return 2 * pow(M_E, x) - 5 * x;
}

double dfunc1(double x) {
    return 2 * pow(M_E, x) - 5;
}

double func2(double x) {
    return pow(M_E, x) + x + 1;
}

double dfunc2(double x) {
    return pow(M_E, x) + 1;
}


int main() {
    // Integration borders
    double a = 1;
    double b = 2;

    // Call trapezoidal rule implementation
    cout << std::setw(10) << std::setprecision(9) << ComputationalMath::TrapezoidIntegration(a, b, func1) << "\n";
    // Call trapezoidal spline method
    cout << std::setw(10) << std::setprecision(9) << ComputationalMath::TrapezoidSplineIntegration(a, b, func1, dfunc1) << "\n";
    // Call Simpson method
    cout << std::setw(10) << std::setprecision(9) << ComputationalMath::SimpsonIntegration(a, b, func1) << "\n";

    return 0;
}
