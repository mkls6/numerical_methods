#include "linear_algebra.hpp"
#include "matrix.hpp"
#include "../include/linear_algebra.hpp"

#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cassert>

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

    for (size_t i = 0; i < matrix.Columns(); i++) {
        for (size_t j = 0; j < matrix.Rows() ; j++) {
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

double LAlgebra::ScalarMultColumnVectors(Matrix& a, Matrix& b){
	if (a.Columns() != 1 || b.Columns() != 1 || a.Rows() != b.Rows())
		throw std::invalid_argument("a.Columns() != 1 or b.Columns() != 1 or a.Rows() != b.Rows()");
	auto res = 0.0;
	for (size_t i = 0; i < a.Rows(); i++)
		res += a[i][0] * b[i][0];
	return res;
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

Matrix LAlgebra::SimpleIterationMethod(Matrix &m, Matrix &b, const double eps, std::ofstream& out){
	auto vec = vector<double>(m.Columns(), 1.0);
	auto x_prev = Matrix(vec);
	auto x_curr = Matrix(1, m.Columns());
	auto stop_criterion = 0.0;
	auto iter = 1;
	auto tau = 0.9 * 2 / LAlgebra::EuclideanNorm(m);
	do {
		auto[x_new, q, residual_norm] = SimpleIterationMethodStep(m, b, x_curr, x_prev, tau);
		auto diff = ~(x_new - x_curr);
		auto tmp = EuclideanNorm(x_prev);
		auto criterion_curr = q * CubicNorm(diff) / (1 - q);
		x_prev = x_curr;
		x_curr = x_new;
		stop_criterion = residual_norm;
		out << "|" << std::setw(5) << iter << "|" << std::setw(14) << tau << "|" << std::setw(14) << q << "|" <<
			std::setw(14) << residual_norm << "|" << std::setw(14) << std::setprecision(7) << std::fixed << criterion_curr << "|" << x_curr;
		iter++;
	} while (stop_criterion > eps);
	return x_curr;
}

tuple<Matrix, double, double> LAlgebra::SimpleIterationMethodStep(Matrix& A, Matrix& b, Matrix& x_curr, Matrix& x_prev, double tau){
	auto answer = Matrix(1, A.Columns());
	for (size_t i = 0; i < A.Rows(); i++) {
		auto sm = 0.0;
		for (size_t j = 0; j < A.Columns(); j++)
			sm += A[i][j] * x_curr[0][j];
		answer[0][i] = x_curr[0][i] + tau * (b[0][i] - sm);
	}
	auto residual = A * ~answer - ~b;
	auto q_numerator = ~(answer - x_curr);
	auto q_denominator = ~(x_curr - x_prev);
	auto q = CubicNorm(q_numerator) / CubicNorm(q_denominator);
	auto residual_norm = EuclideanNorm(residual);
	return { answer, q, residual_norm};
}

Matrix LAlgebra::FastestGradientDescentMethod(Matrix &m, Matrix &b, const double eps, std::ofstream& out){
	auto vec = vector<double>(m.Columns(), 1.0);
	auto x_prev = Matrix(vec);
	auto x_curr = Matrix(1, m.Columns());
	auto stop_criterion = 0.0;
	auto iter = 1;
	auto residual = (m * ~x_curr) - ~b;
	do {
		auto[x_new, residual_curr, tau, q, residual_norm] = FastestGradientDescentMethodStep(m, b, x_curr, x_prev, residual);
		auto diff = ~(x_new - x_curr);
		auto criterion_curr = q * CubicNorm(diff) / (1 - q);
		x_prev = x_curr;
		x_curr = x_new;
		residual = residual_curr;
		stop_criterion = residual_norm;
		out << "|" << std::setw(5) << iter << "|" << std::setw(14) << tau << "|" << std::setw(14) << q << "|" <<
			std::setw(14) << residual_norm << "|" << std::setw(14) << criterion_curr << "|" << x_curr;
		iter++;
	} while (stop_criterion > eps);
	return x_curr;
}

tuple<Matrix, Matrix, double, double, double> LAlgebra::FastestGradientDescentMethodStep(Matrix &A, Matrix &b, Matrix &x_curr, Matrix& x_prev, Matrix& residual_prev){
	auto tmp = A * residual_prev;
	auto tau = ScalarMultColumnVectors(residual_prev, residual_prev) / ScalarMultColumnVectors(tmp, residual_prev);
	auto answer = Matrix(1, A.Columns());
	for (size_t i = 0; i < A.Columns(); i++) {
		auto sm = 0.0;
		for (size_t j = 0; j < A.Columns(); j++)
			sm += A[i][j] * x_curr[0][j];
		answer[0][i] = x_curr[0][i] + tau * (b[0][i] - sm);
	}
	auto residual = ~b - A * (~answer);
	auto q_numerator = ~(answer - x_curr);
	auto q_denominator = ~(x_curr - x_prev);
	auto q = CubicNorm(q_numerator) / CubicNorm(q_denominator);
	auto residual_norm = EuclideanNorm(residual);
	return { answer, residual, tau, q, residual_norm };
}

tuple<double, int> LAlgebra::SorFind(Matrix &m, Matrix &b, double eps, std::ofstream& out){
	auto w = 0.1;
	auto min_iter = -1; 
	auto optim_w = 0.1;
	while (w < 2){
		auto iter = 0;
		auto vec = vector<double>(m.Columns(), 1.0);
		auto x_prev = Matrix(vec);
		auto x_curr = Matrix(1, m.Columns());
		auto residual_norm_criterion = 0.0;
		do{ 
			auto[x_new, q, residual_norm] = SorStep(m, b, x_curr, x_prev, w);
			auto diff = ~(x_new - x_curr);
			auto criterion_curr = q * CubicNorm(diff) / (1 - q);
			x_prev = x_curr;
			x_curr = x_new;
			residual_norm_criterion = residual_norm;
			iter++;
		} while (residual_norm_criterion > eps);
		out << "w= " << w << " itr= " << iter << "\n";
		if (iter < min_iter || min_iter == -1)
		{
			min_iter = iter;
			optim_w = w;
		}
		w += 0.1;
	}
	return { optim_w, min_iter };
}

tuple<Matrix, double, double> LAlgebra::SorStep(Matrix& A, Matrix& b, Matrix& x_curr, Matrix& x_prev, double w)
{
	auto answer = Matrix(1, A.Columns());
	for (size_t i = 0; i < A.Rows(); i++) {
		auto sum_x = 0.0;
		auto sum_x_prev = 0.0;
		for (size_t j = 0; j < i; j++)
			sum_x += A[i][j] * answer[0][j];
		for (size_t j = i + 1; j < A.Columns(); j++)
			sum_x_prev += A[i][j] * x_curr[0][j];
		auto tmp = (b[0][i] - sum_x - sum_x_prev) / A[i][i];
		answer[0][i] = x_curr[0][i] + w * (tmp - x_curr[0][i]);
	}
	auto residual = ~b - A * (~answer);
	auto q_numerator = ~(answer - x_curr);
	auto q_denominator = ~(x_curr - x_prev);
	auto q = CubicNorm(q_numerator) / CubicNorm(q_denominator);
	auto residual_norm = EuclideanNorm(residual);
	return { answer, q, residual_norm };
}

Matrix LAlgebra::Sor(Matrix &m, Matrix &b, const double eps, double w, std::ofstream& out){
	auto vec = vector<double>(m.Columns(), 1.0);
	auto x_prev = Matrix(vec);
	auto x_curr = Matrix(1, m.Columns());
	auto stop_criterion = 0.0;
	auto iter = 1;
	do {
		auto[x_new, q, residual_norm] = SorStep(m, b, x_curr, x_prev, w);
		auto diff = ~(x_new - x_curr);
		auto criterion_curr = q * CubicNorm(diff) / (1 - q);
		x_prev = x_curr;
		x_curr = x_new;
		stop_criterion = residual_norm;
		out << "|" << std::setw(5) << iter << "|" << std::setw(14) << w << "|" << std::setw(14) << q << "|" <<
			std::setw(14) << residual_norm << "|" << std::setw(14) << std::setprecision(7) << std::fixed << criterion_curr << "|" << x_curr;
		iter++;
	} while (stop_criterion > eps);
	return x_curr;
}

Matrix LAlgebra::ConjugateGradientMethod(Matrix &m, Matrix &b, const double eps, std::ofstream& out){
	auto vec = vector<double>(m.Columns(), 1.0);
	auto x_prev = Matrix(vec), x_curr = Matrix(1, m.Columns());
	auto stop_criterion = 0.0;
	auto residual_prev = (m * ~x_curr) - ~b, residual_curr = (m * ~x_curr) - ~b;
	
	auto tmp = m * residual_prev;
	auto tau_prev = ScalarMultColumnVectors(residual_prev, residual_prev) / ScalarMultColumnVectors(tmp, residual_prev);
	
	auto alpha_prev = 0.0;
	auto iter = 1;
	do {
		auto tmp_2 = m * residual_curr;
		auto tau_curr = ScalarMultColumnVectors(residual_curr, residual_curr) / ScalarMultColumnVectors(tmp_2, residual_curr);
		auto alpha_curr = (iter != 1) ? 1 / (1 - tau_curr / tau_prev / alpha_prev * ScalarMultColumnVectors(residual_curr, residual_curr) / ScalarMultColumnVectors(residual_prev, residual_prev)) : 1;

		auto[x_new, residual, q] = ConjugateGradientMethodStep(m, b, x_curr, x_prev, residual_curr, residual_prev, tau_curr, alpha_curr);
		
		residual_prev = residual_curr;
		residual_curr = residual;
		stop_criterion = EuclideanNorm(residual);

		auto diff = ~(x_new - x_curr);
		auto criterion_curr = q * CubicNorm(diff) / (1 - q);
		
		x_prev = x_curr;
		x_curr = x_new;
		tau_prev = tau_curr;
		alpha_prev = alpha_curr;
		
		out << "|" << std::setw(5) << iter << "|" << std::setw(14) << tau_curr << "|" << std::setw(14) << alpha_curr << "|" << std::setw(14) << q << "|" <<
			std::setw(14) << stop_criterion << "|" << std::setw(14) << criterion_curr << "|" << x_curr;
		
		iter++;
	} while (stop_criterion > eps);
	return x_curr;
}

tuple<Matrix, Matrix, double> LAlgebra::ConjugateGradientMethodStep(Matrix &A, Matrix &b, Matrix &x_curr, Matrix& x_prev, Matrix& residual_curr, 
																				Matrix& residual_prev, double tau, double alpha) {
	auto answer = Matrix(1, A.Columns());
	for (size_t i = 0; i < A.Columns(); i++)
		answer[0][i] = alpha * x_curr[0][i] + (1 - alpha) * x_prev[0][i] - tau * alpha * residual_curr[i][0];
	auto residual = A * (~answer) - ~b;
	auto q_numerator = ~(answer - x_curr);
	auto q_denominator = ~(x_curr - x_prev);
	auto q = CubicNorm(q_numerator) / CubicNorm(q_denominator);
	return { answer, residual, q};
}