#ifndef LAB2_LINEAR_ALGEBRA_HPP
#define LAB2_LINEAR_ALGEBRA_HPP

#include <vector>
#include <fstream>
#include "matrix.hpp"

class LAlgebra {
public:
    static double CubicNorm(Matrix &);

    static double OctahedralNorm(Matrix &);

    static double EuclideanNorm(Matrix &);

    static std::vector<double> JacobiEigen(Matrix &);

	static Matrix SimpleIterationMethod(Matrix &, Matrix &, const double, std::ofstream&);

	static Matrix FastestGradientDescentMethod(Matrix &, Matrix &, const double, std::ofstream&);

	static Matrix Sor(Matrix &, Matrix &, const double, double, std::ofstream&);

	static tuple<double, int> SorFind(Matrix &, Matrix &, double, std::ofstream&);

	static Matrix ConjugateGradientMethod(Matrix &, Matrix &, const double, std::ofstream&);

private:
	static tuple<Matrix, double, double> SimpleIterationMethodStep(Matrix &, Matrix &, Matrix &, Matrix &, double);

	static tuple<Matrix, Matrix, double, double, double> FastestGradientDescentMethodStep(Matrix &, Matrix &, Matrix &, Matrix &, Matrix &);
	
	static tuple<Matrix, double, double> SorStep(Matrix &, Matrix &, Matrix &, Matrix &, double);

	static tuple<Matrix, Matrix, double> ConjugateGradientMethodStep(Matrix &, Matrix &, Matrix &, Matrix& , Matrix& , Matrix&, double, double);

	static double ScalarMultColumnVectors(Matrix &, Matrix &);
};

#endif
