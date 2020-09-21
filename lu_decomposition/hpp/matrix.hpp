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
    size_t rows;
    size_t columns;

    Matrix(size_t, size_t);

    Matrix(size_t);

    Matrix(vector<vector<double> > &);

    tuple<Matrix, Matrix, size_t> LUDecompose(optional<vector<int> *>);

    static double LUDeterminant(Matrix &,
                                Matrix &,
                                size_t);

    static double LUSolve(Matrix &,
                          Matrix &);

    static Matrix *ReadMatrix(istream &);

    vector<double> operator[](const size_t &index);

    friend ostream &operator<<(ostream &, const Matrix &);

private:
    vector<vector<double> > data;
};


#endif
