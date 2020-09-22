#include "../hpp/linear_algebra.hpp"
#include <cmath>

double LAlgebra::CubicNorm(Matrix &matrix) {
    double max_sum = 0;
    double tmp_sum = 0;

    for (size_t i = 0; i < matrix.Rows(); i++) {
        for (size_t j = 0; j < matrix.Columns(); j++) {
            tmp_sum += fabs(matrix[i][j]);
        }

        if (tmp_sum > max_sum)
            max_sum = tmp_sum;
        tmp_sum = 0;
    }

    return max_sum;
}

double LAlgebra::OctahedralNorm(Matrix &matrix) {
    double max_sum = 0;
    double tmp_sum = 0;

    for (size_t i = 0; i < matrix.Rows(); i++) {
        for (size_t j = 0; j < matrix.Columns(); j++) {
            tmp_sum += fabs(matrix[j][i]);
        }
        if (tmp_sum > max_sum)
            max_sum = tmp_sum;
        tmp_sum = 0;
    }

    return max_sum;
}

double LAlgebra::EuclideanNorm(Matrix &) {
    return 0;
}

double LAlgebra::JacobiRotation() {
    return 0;
}
