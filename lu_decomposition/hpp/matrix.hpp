#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <memory>
#include <utility>
#include <optional>
#include <vector>
#include <ostream>
#include <istream>

using std::tuple, std::optional, std::vector, std::ostream, std::istream;

class Matrix {
public:
    Matrix(size_t, size_t);

    Matrix(size_t);

    Matrix(vector<vector<double> > &);

    tuple<Matrix, Matrix, Matrix, size_t> LUDecompose();

    static double LUDeterminant(Matrix &,
                                Matrix &,
                                size_t);

    static vector<double> LUSolve(Matrix &,
                                  Matrix &,
                                  vector<double> &);

    static Matrix LUInverseMatrix(Matrix &, Matrix &, Matrix &);

    static Matrix *ReadMatrix(istream &);

    Matrix operator+(const Matrix &);

    Matrix operator-(const Matrix &);

    Matrix operator*(const Matrix &);

    Matrix operator~();

    vector<double> &operator[](const size_t &index);

    friend ostream &operator<<(ostream &, const Matrix &);

private:
    size_t rows;
    size_t columns;

    vector<vector<double> > data;
};


#endif
