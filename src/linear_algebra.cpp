#include "../include/linear_algebra.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>

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

double LAlgebra::EuclideanNorm(Matrix &m) {
    auto mT = ~m;
    auto mul = mT * m;
    auto eigenValues = JacobiEigen(mul);

    return sqrt(*(std::max_element(eigenValues.begin(), eigenValues.end())));
}

double offValue(Matrix &m) {
    double sum = 0;

    for (size_t i = 0; i < m.Rows(); i++) {
        for (size_t j = i + 1; j < m.Columns(); j++) {
            sum += m[i][j] * m[i][j] + m[j][i] * m[j][i];
        }
    }

    return sqrt(sum);
}

std::vector<double> LAlgebra::JacobiEigen(Matrix &m) {
    double threshold = 1e-10;
    double pivot;
    double currentOffValue, prevOffValue;
    std::pair<size_t, size_t> pivotPos;
    auto rotation = Matrix(m.Rows());
    auto copy = Matrix(m);

    for (size_t i = 0; i < rotation.Rows(); i++)
        rotation[i][i] = 1;

    prevOffValue = currentOffValue = offValue(copy);
    prevOffValue += 10;

    while (currentOffValue > threshold && fabs(currentOffValue - prevOffValue) > threshold) {
        for (size_t i = 0; i < copy.Rows(); i++) {
            // Reset the pivot
            pivot = 0;
            pivotPos = {i, i};

            for (size_t j = 0; j < copy.Columns(); j++) {
                // Skip diagonal items
                if (i == j)
                    continue;

                // Find max off-diagonal item in current row
                if (fabs(copy[i][j]) > fabs(pivot)) {
                    pivot = copy[i][j];
                    pivotPos = {i, j};
                }
            }

            // Calculate cos(theta) and s(theta)
            auto[p, q] = pivotPos;
            double angleTmp = 2 * copy[p][q] / (copy[p][p] - copy[q][q]);
            auto sign = (angleTmp > 0) - (angleTmp < 0);
            double c = sign * sqrt(0.5 * (1 + 1 / sqrt(1 + angleTmp * angleTmp)));
            double s = sign * sqrt(0.5 * (1 - 1 / sqrt(1 + angleTmp * angleTmp)));

            // Set up Jacobi rotation matrix
            rotation[p][p] = rotation[q][q] = c;
            rotation[p][q] = s;
            rotation[q][p] = -rotation[p][q];

            // Perform rotation
            copy = ~rotation * copy * rotation;

            // Reset rotation matrix
            rotation[p][q] = rotation[q][p] = 0;
            rotation[p][p] = rotation[q][q] = 1;
        }

        prevOffValue = currentOffValue;
        currentOffValue = offValue(copy);
    }

    vector<double> eigenValues(copy.Rows());

    for (size_t i = 0; i < copy.Rows(); i++) {
        eigenValues[i] = copy[i][i];
    }

    return eigenValues;
}

Matrix LAlgebra::NLSimpleIterSolve(double x0,
                                   double y0,
                                   vector<std::function<double(double, double)>> &functions,
                                   vector<std::function<double(double)>> &stdFunctions,
                                   vector<vector<std::function<double(double, double)>>> &derivatives) {
    double x = x0;
    double y = y0;
    double residualNorm = 0;
    double q = 0;

    auto f1 = functions[0];
    auto f2 = functions[1];
    auto f1n = stdFunctions[0];
    auto f2n = stdFunctions[1];

    size_t iter = 0;

    Matrix Jacobian(2, 2);
    for (size_t i = 0; i < Jacobian.Rows(); i++) {
        for (size_t j = 0; j < Jacobian.Columns(); j++) {
            Jacobian[i][j] = derivatives[i][j](x0, y0);
        }
    }

    std::cout << "Jacobian:\n" << Jacobian << "\n";
    std::cout << "Jacobian norm:\n" << CubicNorm(Jacobian) << "\n";
    std::cout << std::setw(18) << "Itr"
              << std::setw(18)  << "x"
              << std::setw(18) << "y"
              << std::setw(18) << "Residual norm"
              << std::setw(18) << "F1"
              << std::setw(18) << "F2"
              << std::setw(18) << "Jacobian norm"
              << "\n";

    do {
        iter++;

        x = f1n(y0);
        y = f2n(x0);

        for (size_t i = 0; i < Jacobian.Rows(); i++) {
            for (size_t j = 0; j < Jacobian.Columns(); j++) {
                Jacobian[i][j] = derivatives[i][j](x, y);
            }
        }
        q = CubicNorm(Jacobian);
        residualNorm = std::max(fabs(x - x0), fabs(y - y0)) * (1 - q) / q;
        std::cout << std::setprecision(10)
                  << std::setw(18) << iter
                  << std::setw(18) << x
                  << std::setw(18) << y
                  << std::setw(18) << residualNorm
                  << std::setw(18) << f1(x, y)
                  << std::setw(18) << f2(x, y)
                  << std::setw(18) << q
                  << "\n";

        x0 = x;
        y0 = y;
    } while (residualNorm > eps);

    vector<double> result({x, y});
    return result;
}
