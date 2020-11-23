#include "../include/computational_math.hpp"

#include <cmath>
using std::fabs;

double ComputationalMath::TrapezoidIntegration(double left,
                                               double right,
                                               const function<double(double)>& f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width = 1;
    // Computed integral value
    double previousValue = 0;
    double currentValue = 0;
    // Runge rule coefficient for trapezoid rule
    double theta = 1 / 3.;

    do {
        double x1 = left;
        previousValue = currentValue;
        currentValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x2 = x1 + width;
            currentValue += 0.5 * (x2 - x1) * (f(x1) + f(x2));
            x1 = x2;
        }

        segmentsNumber *= 2;
    } while (theta * fabs(currentValue - previousValue) > EPS);

    return currentValue;
}

double ComputationalMath::TrapezoidSplineIntegration(double left,
                                                     double right,
                                                     const function<double(double)>& f,
                                                     const function<double(double)>& df) {
    return 0;
}

double ComputationalMath::SimpsonIntegration(double left,
                                             double right,
                                             const function<double(double)> &f,
                                             const function<double(double)> &df) {
    return 0;
}

double ComputationalMath::Newton3Integration(double, double, const function<double(double)>&,
                                             const function<double(double)>&) {
    return 0;
}
