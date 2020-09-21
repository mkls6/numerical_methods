#include <cassert>
#include <string>
#include <iomanip>
#include <cmath>
#include <iostream>

#include "../hpp/matrix.hpp"

ostream &operator<<(ostream &out, const Matrix &matrix) {
    out << "Matrix:\n" << std::fixed;
    for (size_t r = 0; r < matrix.rows; r++) {
        for (size_t c = 0; c < matrix.columns; c++) {
            out << std::setw(14) << std::setprecision(8) << matrix.data[r][c] << " ";
        }
        out << '\n';
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

double Matrix::LUSolve(Matrix &L, Matrix &U) {
    return 0;
}

size_t find_max(vector<vector<double> > &src, size_t i, size_t j) {
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

tuple<Matrix, Matrix, size_t> Matrix::LUDecompose(optional<vector<int> *> permutations) {
    assert(this->rows == this->columns && "Unable to perform LU decomposition on non-square matrix");

    size_t n = this->rows;
    size_t permutationsCount = 0;

    vector<vector<double> > L(n, vector<double>(n, 0));
    vector<vector<double> > U(this->data);

    for (size_t i = 0; i < n; i++) {
        std::cout << "Step " << i + 1 << "\n";
        auto max_pos = find_max(U, i, i);
        // Check if max value is zero or really close to zero
        if (fabs(U[i][max_pos]) < 1e-8)
            throw std::invalid_argument("Max element is within [0:1e-8]. Aborting.");

        // Swap lines if needed
        if (max_pos != i) {
            std::swap(L[i], L[max_pos]);
            std::swap(U[i], U[max_pos]);
            permutationsCount++;

            if (permutations.has_value())
                std::swap(permutations.value()[i], permutations.value()[max_pos]);
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

        auto temp = U[i][i];
        for (size_t k = i; k < n; k++) {
            U[i][k] /= temp;
        }

        std::cout << "L:\n" << L << "\nU:\n" << U << "\n";
    }

    return {U, L, permutationsCount};
}

Matrix::Matrix(size_t rows, size_t columns) {
    this->rows = rows;
    this->columns = columns;

    this->data.resize(rows);
    for (auto i : this->data) {
        i.resize(columns, 0);
    }
}

Matrix::Matrix(size_t n) : Matrix(n, n) {}

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

vector<double> Matrix::operator[](const size_t &index) {
    return this->data[index];
}
