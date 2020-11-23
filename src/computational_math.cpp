#include "../include/computational_math.hpp"

#include <cmath>
using std::fabs, std::sqrt;

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
        width = (right - left) / segmentsNumber;

        previousValue = currentValue;
        currentValue = width * 0.5 * (f(left) + f(right)) + pow(width, 2) / 12. * (df(left) - df(right));
        double innerSum = 0;

        double x = left;
        for (int step = 1; step < segmentsNumber; step++) {
            x += width;
            innerSum += f(x);
        }
        currentValue += width * innerSum;
        segmentsNumber *= 2;

    } while (theta * fabs(currentValue - previousValue) > EPS);

    return currentValue;
}

double ComputationalMath::SimpsonIntegration(double left,
                                             double right,
                                             const function<double(double)> &f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width = 1;
    // Computed integral value
    double previousValue = 0;
    double currentValue = 0;
    // Runge rule coefficient for Simpson method
    double theta = 1 / 15.;

    do {
        double x1 = left;
        previousValue = currentValue;
        currentValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x2 = x1 + width;
            currentValue += (x2 - x1) / 6. * (f(x1) + 4.0 * f(0.5 * (x1 + x2)) + f(x2));
            x1 = x2;
        }

        segmentsNumber *= 2;
    } while (theta * fabs(currentValue - previousValue) > EPS);

    return currentValue;
}

double ComputationalMath::Gauss3Integration(double left,
                                            double right,
                                            const function<double(double)>& f) {
    // Initial values
    double segmentsNumber = 1;
    // Width of the segment
    double width = 1;
    // Computed integral value
    double previousValue = 0;
    double currentValue = 0;
    // Runge rule coefficient for Newton method
    double theta = 1 / 63.;

    do {
        previousValue = currentValue;
        currentValue = 0;
        width = (right - left) / segmentsNumber;

        for (int step = 0; step < segmentsNumber; step++) {
            double x = left + 0.5 * (1 + 2 * step) * width;
            currentValue += width / 2. * (f(x - width * sqrt(3 / 5.) / 2.) * (5 / 9.) + (8 / 9.)
                    * f(x) + (5 / 9.) * f(x + width * sqrt(3 / 5.) / 2.));
        }

        segmentsNumber *= 2;
    } while (theta * fabs(currentValue - previousValue) > EPS);

    return currentValue;
}
