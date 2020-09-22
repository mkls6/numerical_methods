#ifndef LAB2_LINEAR_ALGEBRA_HPP
#define LAB2_LINEAR_ALGEBRA_HPP

#include <vector>

#include "matrix.hpp"

class LAlgebra {
public:
    static double CubicNorm(Matrix &);

    static double OctahedralNorm(Matrix &);

    static double EuclideanNorm(Matrix &);

    static std::vector<double> JacobiEigen(Matrix &);
};

#endif
