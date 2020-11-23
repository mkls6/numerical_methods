#ifndef COMPUTATIONAL_MATH_HPP
#define COMPUTATIONAL_MATH_HPP


#include <functional>

using std::function;

class ComputationalMath {
    static double TrapezoidIntegration(double,
                                       double,
                                       const function<double(double, double)>&,
                                       const function<double(double, double)>&);

    static double TrapezoidSplineIntegration(double,
                                             double,
                                             const function<double(double, double)>&,
                                             const function<double(double, double)>&);

    static double SimpsonIntegration(double,
                                     double,
                                     const function<double(double, double)>&,
                                     const function<double(double, double)>&);

    static double Newton3Integration(double,
                                     double,
                                     const function<double(double, double)>&,
                                     const function<double(double, double)>&);
};


#endif
