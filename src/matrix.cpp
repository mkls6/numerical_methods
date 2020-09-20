#include "../hpp/matrix.hpp"
#include <cassert>

ostream &operator<<(ostream &out, const Matrix &matrix) {
    out << "Matrix:\n";
    for (size_t r = 0; r < matrix.rows; r++) {
        for (size_t c = 0; c < matrix.columns; c++) {
            out << matrix.data[r][c] << " ";
        }
        out << '\n';
    }

    return out;
}

double Matrix::LUDeterminant(Matrix *, Matrix *) {
    return 0;
}

double Matrix::LUSolve(Matrix *, Matrix *) {
    return 0;
}

size_t find_max(vector<vector<double> > &src, size_t i, size_t j) {
    pair<size_t, size_t> max_pos = {i, j};
    double max_element = src[i][j];

    for (size_t k = j + 1; k < src[0].size(); k++) {
        if (src[i][k] > max_element) {
            max_element = src[i][k];
            max_pos = {i, k};
        }
    }

    return max_pos.second;
}

pair<Matrix, Matrix> Matrix::LUDecompose(optional<vector<int>> &permutations) {
    assert(this->rows == this->columns && "Unable to perform LU decomposition on non-square matrix");

    size_t n = this->rows;

    vector<vector<double> > L(n, vector<double>(n, 0));
    vector<vector<double> > U(n, vector<double>(n, 0));
    vector<vector<double> > originalCopy(n, vector<double>(n, 0));
    std::copy(this->data.begin(), this->data.end(), originalCopy);

    for (size_t i = 0; i < n; i++) {
        std::swap(originalCopy[find_max(originalCopy, i, i)], originalCopy[i]);
        if (originalCopy[i][i] == 0)
            throw std::invalid_argument("Invalid matrix");

        for (size_t j = 0; j <= i; j++) {
            L[i][j] = this->data[i][j];

            for (size_t k = 1; long(k) < long(j) - 1; k++) {
                L[i][j] -= U[k][i] * L[i][k];
            }
        }
        for (size_t j = i + 1; j < n; j++) {
            U[i][j] = this->data[i][j];

            for (size_t k = 1; long(k) < long(i) - 1; k++) {
                U[i][j] -= U[k][j] * L[i][k];
            }

            U[i][j] /= L[i][i];
        }
    }

    return {U, L};
}

Matrix::Matrix(size_t rows, size_t columns) {
    this->rows = rows;
    this->columns = columns;

    this->data.resize(rows);
    for (auto i : this->data) {
        i.resize(columns, 0);
    }
}

Matrix::Matrix(size_t n) : Matrix(n, n) { }

Matrix::Matrix(vector<vector<double>> &data) {
    assert(!data.empty() && !data[0].empty() && "Empty vectors were provided.");

    this->data = data;
    this->rows = data.size();
    this->columns = data[0].size();
}
