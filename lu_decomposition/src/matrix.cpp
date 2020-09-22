#include <cassert>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <iostream>

#include "../hpp/matrix.hpp"

ostream &operator<<(ostream &out, const Matrix &matrix) {
    for (size_t r = 0; r < matrix.rows; r++) {
        for (size_t c = 0; c < matrix.columns; c++) {
            out << std::setw(16) << matrix.data[r][c];
        }
        out << "\n";
    }

    return out;
}

// Actually, we need only one of the given matrices
// as one of them always has det = 1. But we'll calculate
// both in case of different LU decomposition method's been used
double Matrix::LUDeterminant(Matrix &L, Matrix &U, size_t permutationCount) {
    assert(L.rows == L.columns && "Non-square L given.");
    assert(U.rows == U.columns && "Non-square U given.");
    assert(L.rows == U.rows && "L and U have different sizes.");

    double detL = 1, detU = 1;

    for (size_t i = 0; i < L.rows; i++) {
        detL *= L[i][i];
        detU *= U[i][i];
    }

    return detL * detU * (permutationCount % 2 ? -1 : 1);
}

vector<double> Matrix::LUSolve(Matrix &L, Matrix &U, vector<double> &b) {
    assert(L.rows == L.columns && "Non-square L given.");
    assert(U.rows == U.columns && "Non-square U given.");
    assert(L.rows == U.rows && "L and U have different sizes.");

    vector<double> X(L.rows, 0);
    vector<double> Y(L.rows, 0);

    for (size_t i = 0; i < L.rows; i++) {
        double tmp = b[i];
        for (size_t j = 0; j < i; j++)
            tmp -= L[i][j] * Y[j];
        Y[i] = tmp / L[i][i];
    }
    for (long i = long(L.rows) - 1; i >= 0; i--) {
        double tmp = Y[i];
        for (size_t j = i + 1; j < L.rows; j++)
            tmp -= U[i][j] * X[j];
        X[i] = tmp / U[i][i];
    }

    return X;
}

size_t findMax(vector<vector<double> > &src, size_t i, size_t j) {
    std::pair<size_t, size_t> max_pos = {i, j};
    double max_element = src[i][j];

    for (size_t k = i + 1; k < src[0].size(); k++) {
        if (fabs(src[k][j]) > fabs(max_element)) {
            max_element = src[k][j];
            max_pos = {k, j};
        }
    }

    return max_pos.first;
}

// Return L, U, P, number of permutations
tuple<Matrix, Matrix, Matrix, size_t> Matrix::LUDecompose() {
    assert(this->rows == this->columns && "Unable to perform LU decomposition on non-square matrix");

    size_t n = this->rows;
    size_t permutationsCount = 0;

    vector<vector<double> > L(n, vector<double>(n, 0));
    vector<vector<double> > U(this->data);
    vector<vector<double> > P(n, vector<double>(n, 0));

    for (size_t i = 0; i < n; i++)
        P[i][i] = 1;

    for (size_t i = 0; i < n; i++) {
        std::cout << "Step " << i + 1 << "\n";
        auto max_pos = findMax(U, i, i);
        // Check if max value is zero or really close to zero
        if (fabs(U[i][max_pos]) < 1e-8)
            throw std::invalid_argument("Max element is within [0:1e-8]. Aborting.");

        // Swap lines if needed
        if (max_pos != i) {
            std::cout << "Swapping lines " << i + 1 << " and " << max_pos + 1 << "\n";
            std::swap(L[i], L[max_pos]);
            std::swap(U[i], U[max_pos]);
            std::swap(P[i], P[max_pos]);
            permutationsCount++;
        }

        // L[i][j] = U[i][j] - SUM(L[i][k] * U[k][i]), k = 1..(j-1)
        for (size_t j = i; j < n; j++) {
            L[j][i] = U[j][i];
        }

        // Perform Gaussian elimination
        for (size_t k = i + 1; k < n; k++) {
            auto delta = U[k][i] / U[i][i];
            for (size_t j = i; j < n; j++) {
                U[k][j] -= U[i][j] * delta;
            }
        }

        // Divide line by pivot value
        auto temp = U[i][i];
        for (size_t k = i; k < n; k++) {
            U[i][k] /= temp;
        }

        std::cout << "L:\n" << L << "\nU:\n" << U << "\n";
    }

    return {L, U, P, permutationsCount};
}

Matrix::Matrix(size_t rows, size_t columns) {
    this->rows = rows;
    this->columns = columns;

    this->data.resize(rows);
    for (auto &i : this->data) {
        i.resize(columns, 0);
    }
}

Matrix::Matrix(size_t n) : Matrix(n, n) {}

Matrix::Matrix(vector<double> &v) {
    this->rows = 1;
    this->columns = v.size();

    this->data.push_back(v);
}

Matrix::Matrix(vector<vector<double>> &data) {
    assert(!data.empty() && !data[0].empty() && "Empty vectors were provided.");

    this->data = data;
    this->rows = data.size();
    this->columns = data[0].size();
}

Matrix *Matrix::ReadMatrix(istream &in) {
    vector<vector<double> > lines;
    vector<double> line;
    std::string str_line;

    while (std::getline(in, str_line)) {
        std::istringstream iss(str_line);
        std::string tmp;

        while (iss >> tmp) {
            line.push_back(std::stod(tmp));
        }
        lines.push_back(line);
        line.clear();
    }

    return new Matrix(lines);
}

Matrix *Matrix::ReadMatrix(istream &in, size_t n) {
    return Matrix::ReadMatrix(in, n, n);
}

Matrix *Matrix::ReadMatrix(istream &in, size_t rows, size_t columns) {
    auto *m = new Matrix(rows, columns);

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < columns; j++) {
            in >> (*m)[i][j];
        }
    }

    return m;
}


vector<double> &Matrix::operator[](const size_t &index) {
    return this->data[index];
}

Matrix Matrix::operator+(const Matrix &m) {
    if (this->rows != m.rows || this->columns != m.columns)
        throw std::invalid_argument("Unequal matrix sizes");

    auto newMatrix = Matrix(m.rows, m.columns);

    for (size_t i = 0; i < m.rows; i++) {
        for (size_t j = 0; j < m.columns; j++)
            newMatrix[i][j] = this->data[i][j] + m.data[i][j];
    }

    return newMatrix;
}

Matrix Matrix::operator-(const Matrix &m) {
    if (this->rows != m.rows || this->columns != m.columns)
        throw std::invalid_argument("Unequal matrix sizes");

    auto newMatrix = Matrix(m.rows, m.columns);

    for (size_t i = 0; i < m.rows; i++) {
        for (size_t j = 0; j < m.columns; j++)
            newMatrix[i][j] = this->data[i][j] - m.data[i][j];
    }

    return newMatrix;
}

Matrix Matrix::operator*(const Matrix &m) {
    if (this->columns != m.rows)
        throw std::invalid_argument("A x B: A column count must be equal to B row count.");

    auto newMatrix = Matrix(this->rows, m.columns);

    for (size_t i = 0; i < this->rows; i++) {
        for (size_t j = 0; j < m.columns; j++) {
            for (size_t k = 0; k < this->columns; k++) {
                newMatrix[i][j] += this->data[i][k] * m.data[k][j];
            }
        }
    }

    return newMatrix;
}

Matrix Matrix::LUInverseMatrix(Matrix &L, Matrix &U, Matrix &P) {
    assert(L.rows == U.rows && L.columns == U.columns && L.rows == L.columns);

    auto inverseMatrix = Matrix(L.rows, L.columns);

    for (size_t i = 0; i < L.rows; i++) {
        inverseMatrix[i] = Matrix::LUSolve(L, U, P[i]);
    }

    return ~inverseMatrix;
}

Matrix Matrix::operator~() {
    auto newMatrix = Matrix(this->columns, this->rows);
    for (size_t i = 0; i < this->rows; i++) {
        for (size_t j = 0; j < this->columns; j++) {
            newMatrix[j][i] = this->data[i][j];
        }
    }
    return newMatrix;
}

size_t Matrix::Rows() const {
    return this->rows;
}

size_t Matrix::Columns() const {
    return this->columns;
}
