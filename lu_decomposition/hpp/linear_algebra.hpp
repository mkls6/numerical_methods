#ifndef LAB2_LINEAR_ALGEBRA_HPP
#define LAB2_LINEAR_ALGEBRA_HPP

#include "matrix.hpp"

class LAlgebra {
public:
    static double CubicNorm(Matrix &);

    static double OctahedralNorm(Matrix &);

    static double EuclideanNorm(Matrix &);

    static double JacobiRotation();
};

#endif
