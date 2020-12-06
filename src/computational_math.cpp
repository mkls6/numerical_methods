#include "../include/computational_math.hpp"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

using std::fabs, std::sqrt, std::log, std::cout, std::to_string;

const int decimals = 10;

double KValue(double prev, double current, double next) {
    return log((next - prev) / (current - prev) - 1) / log(0.5);
}

void PrintHeader() {
    cout << std::setw(6) << "N"
         << std::fixed
         << std::setw(decimals + 2) << "h"
         << std::setw(decimals + 2) << std::setprecision(decimals) << "Integral"
         << std::setw(decimals + 2) << "Error"
         << std::setw(decimals + 2) << "k\n";
}

double ComputationalMath::TrapezoidIntegration(double left,
                                               double right,
                                               const function<double(double)> &f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width;
    // Computed integral value
    double previousValue;
    double currentValue = 0;
    double nextValue = 0;
    // Runge rule coefficient for trapezoid rule
    double theta = 1 / 3.;
    double error;
    double leftF = f(left);

    // Print output header
    PrintHeader();

    do {
        double x1 = left;
        double previousFuncValue = leftF;
        previousValue = currentValue;
        currentValue = nextValue;

        nextValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x2 = x1 + width;
            double funcValue = f(x2);
            nextValue += 0.5 * (x2 - x1) * (previousFuncValue + funcValue);
            previousFuncValue = funcValue;
            x1 = x2;
        }

        error = theta * fabs(nextValue - currentValue);

        cout << std::setprecision(0) << std::setw(6) << segmentsNumber
             << std::setprecision(decimals - 2)
             << std::fixed << std::setw(decimals + 2) << width
             << std::setw(decimals + 2) << std::fixed << std::setprecision(decimals - 2) << nextValue
             << std::scientific << std::setprecision(decimals - 6) << std::setw(decimals + 2) << error
             << std::fixed << std::setprecision(decimals - 1) << std::setw(decimals + 2)
             << (segmentsNumber >= 4 ? to_string(KValue(previousValue, currentValue, nextValue)) : "")
             << "\n";

        segmentsNumber *= 2;
    } while (error > EPS);

    return currentValue;
}

double ComputationalMath::TrapezoidSplineIntegration(double left,
                                                     double right,
                                                     const function<double(double)> &f,
                                                     const function<double(double)> &df) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width;
    // Computed integral value
    double previousValue;
    double currentValue = 0;
    double nextValue = 0;
    // Runge rule coefficient for trapezoid rule
    double theta = 1 / 3.;
    double error;

    double leftF = f(left);
    double rightF = f(right);
    double leftDF = df(left);
    double rightDF = df(right);

    do {
        width = (right - left) / segmentsNumber;

        previousValue = currentValue;
        currentValue = nextValue;

        nextValue = width * 0.5 * (leftF + rightF) + pow(width, 2) / 12. * (leftDF - rightDF);
        double innerSum = 0;

        double x = left;
        for (int step = 1; step < segmentsNumber; step++) {
            x += width;
            innerSum += f(x);
        }
        nextValue += width * innerSum;

        error = theta * fabs(nextValue - currentValue);

        cout << std::setprecision(0) << std::setw(6) << segmentsNumber
             << std::setprecision(decimals - 2)
             << std::fixed << std::setw(decimals + 2) << width
             << std::setw(decimals + 2) << std::fixed << std::setprecision(decimals - 2) << nextValue
             << std::scientific << std::setprecision(decimals - 6) << std::setw(decimals + 2) << error
             << std::fixed << std::setprecision(decimals - 1) << std::setw(decimals + 2)
             << (segmentsNumber >= 4 ? to_string(KValue(previousValue, currentValue, nextValue)) : "")
             << "\n";

        segmentsNumber *= 2;
    } while (error > EPS);

    return currentValue;
}

double ComputationalMath::SimpsonIntegration(double left,
                                             double right,
                                             const function<double(double)> &f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width;
    // Computed integral value
    double previousValue;
    double currentValue = 0;
    double nextValue = 0;
    // Runge rule coefficient for Simpson method
    double theta = 1 / 15.;
    double error;
    double leftF = f(left);

    do {
        double x1 = left;
        previousValue = currentValue;
        currentValue = nextValue;
        double previousFuncValue = leftF;

        nextValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x2 = x1 + width;
            double newFuncValue = f(x2);
            nextValue += (x2 - x1) / 6. * (previousFuncValue + 4.0 * f(0.5 * (x1 + x2)) + newFuncValue);
            previousFuncValue = newFuncValue;
            x1 = x2;
        }

        error = theta * fabs(nextValue - currentValue);

        cout << std::setprecision(0) << std::setw(6) << segmentsNumber
             << std::setprecision(decimals - 2)
             << std::fixed << std::setw(decimals + 2) << width
             << std::setw(decimals + 2) << std::fixed << std::setprecision(decimals - 2) << nextValue
             << std::scientific << std::setprecision(decimals - 6) << std::setw(decimals + 2) << error
             << std::fixed << std::setprecision(decimals - 1) << std::setw(decimals + 2)
             << (segmentsNumber >= 4 ? to_string(KValue(previousValue, currentValue, nextValue)) : "")
             << "\n";

        segmentsNumber *= 2;
    } while (error > EPS);

    return currentValue;
}

double ComputationalMath::GaussianQuadratureIntegration(double left,
                                                        double right,
                                                        const function<double(double)> &f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width;
    // Computed integral value
    double previousValue;
    double currentValue = 0;
    double nextValue = 0;
    // Runge rule coefficient for Newton method
    double theta = 1 / 63.;
    double error;

    do {
        previousValue = currentValue;
        currentValue = nextValue;

        nextValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x = left + 0.5 * (1 + 2 * step) * width;
            nextValue += width / 2. * (f(x - width * sqrt(3 / 5.) / 2.) * (5 / 9.) + (8 / 9.)
                                                                                     * f(x) +
                                       (5 / 9.) * f(x + width * sqrt(3 / 5.) / 2.));
        }

        error = theta * fabs(nextValue - currentValue);

        cout << std::setprecision(0) << std::setw(6) << segmentsNumber
             << std::setprecision(decimals - 2)
             << std::fixed << std::setw(decimals + 2) << width
             << std::setw(decimals + 2) << std::fixed << std::setprecision(decimals - 2) << nextValue
             << std::scientific << std::setprecision(decimals - 6) << std::setw(decimals + 2) << error
             << std::fixed << std::setprecision(decimals - 1) << std::setw(decimals + 2)
             << (segmentsNumber >= 4 ? to_string(KValue(previousValue, currentValue, nextValue)) : "")
             << "\n";

        segmentsNumber *= 2;
    } while (error > EPS);

    return currentValue;
}
