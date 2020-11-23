#ifndef COMPUTATIONAL_MATH_HPP
#define COMPUTATIONAL_MATH_HPP


#include <functional>

using std::function;

class ComputationalMath {
public:
    static constexpr double EPS = 1e-8;
    static double TrapezoidIntegration(double,
                                       double,
                                       const function<double(double)>&);

    static double TrapezoidSplineIntegration(double,
                                             double,
                                             const function<double(double)>&,
                                             const function<double(double)>&);

    static double SimpsonIntegration(double,
                                     double,
                                     const function<double(double)>&);

    static double Gauss3Integration(double,
                                    double,
                                    const function<double(double)>&);
};


#endif
