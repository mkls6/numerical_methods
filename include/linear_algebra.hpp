#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <vector>
#include <functional>

#include "matrix.hpp"

class LAlgebra {
public:
    static constexpr double eps = 1e-4;
    static double CubicNorm(Matrix &);

    static double OctahedralNorm(Matrix &);

    static double EuclideanNorm(Matrix &);

    static std::vector<double> JacobiEigen(Matrix &);

    static Matrix NLSimpleIterSolve(double,
                                    double,
                                    vector<std::function<double(double, double)> > &,
                                    vector<std::function<double(double)>> &,
                                    vector<vector<std::function<double(double, double)> > > &);

    static Matrix NLNewtonSolve(double,
                                double,
                                vector<std::function<double(double, double)>> &,
                                vector<vector<std::function<double(double, double)> > > &,
                                vector<vector<std::function<double(double, double)> > > &);

    static Matrix NLGradientDescentSolve(double,
                                         double,
                                         double,
                                         double,
                                         vector<std::function<double(double, double)>> &,
                                         vector<std::function<double(double, double)>> &,
                                         vector<vector<std::function<double(double, double)>>> &);
};

#endif
