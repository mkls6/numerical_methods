#include <cmath>
#include <iostream>
#include "../include/computational_math.hpp"

using std::pow, std::cout;

double func1(double x) {
    return 2 * pow(M_E, x) - 5 * x;
}

double func2(double x) {
    return pow(M_E, x) + x + 1;
}


int main() {
    // Call trapezoidal rule implementation
    cout << ComputationalMath::TrapezoidIntegration(1, 2, func1) << "\n";

    return 0;
}
