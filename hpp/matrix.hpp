#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <memory>
#include <utility>
#include <optional>
#include <vector>
#include <ostream>

using std::array,
std::pair,
std::optional,
std::vector,
std::ostream;

class Matrix {
public:
    size_t rows;
    size_t columns;

    Matrix(size_t, size_t);
    Matrix(size_t);
    Matrix(vector<vector<double> > &);

    pair<Matrix, Matrix> LUDecompose(optional<vector<int> > &);
    static double LUDeterminant(Matrix *,
                                Matrix *);
    static double LUSolve(Matrix *,
                          Matrix *);

    friend ostream &operator<<(ostream &, const Matrix &);

private:
    vector<vector<double> > data;
};


#endif
