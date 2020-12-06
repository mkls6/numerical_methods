#include <cmath>
#include <iostream>
#include <iomanip>
#include "../include/computational_math.hpp"

using std::pow, std::cout;

int numberOfFunctionCalls = 0;

double func1(double x) {
    numberOfFunctionCalls++;
    return 2 * pow(M_E, x) - 5 * x;
}

double df1(double x) {
    return 2 * pow(M_E, x) - 5;
}

double func2(double x) {
    numberOfFunctionCalls++;
    return pow(M_E, x) + x + 1;
}

double df2(double x) {
    return pow(M_E, x) + 1;
}


int main() {
    // Integration borders
    double a = 1;
    double b = 2;

    // Call trapezoidal rule implementation
    cout << "\nTrapezoid\n";
    auto value = ComputationalMath::TrapezoidIntegration(a, b, func2);
    cout << "Result: " << std::setprecision(16) << std::defaultfloat << value << "\nFunction calls: "
         << numberOfFunctionCalls << "\n";
    numberOfFunctionCalls = 0;

    // Call trapezoidal spline method
    cout << "\nTrapezoid with spline\n";
    value = ComputationalMath::TrapezoidSplineIntegration(a, b, func2, df2);
    cout << "Result: " << std::setprecision(16) << std::defaultfloat << value << "\nFunction calls: "
         << numberOfFunctionCalls << "\n";
    numberOfFunctionCalls = 0;

    // Call Simpson method
    cout << "\nSimpson\n";
    value = ComputationalMath::SimpsonIntegration(a, b, func2);
    cout << "Result: " << std::setprecision(16) << std::defaultfloat << value << "\nFunction calls: "
         << numberOfFunctionCalls << "\n";
    numberOfFunctionCalls = 0;

    // Call Newton3 method
    cout << "\nGaussian quadrature\n";
    value = ComputationalMath::GaussianQuadratureIntegration(a, b, func2);
    cout << "Result: " << std::setprecision(16) << std::defaultfloat << value << "\nFunction calls: "
         << numberOfFunctionCalls << "\n";
    numberOfFunctionCalls = 0;

    return 0;
}
