#include <cassert>
#include "../hpp/linear_algebra.hpp"

size_t find_max(vector<vector<double> > &matrix,
                const size_t start_row,
                const size_t start_column) {
    assert(start_row < matrix.size() &&
           start_column < matrix.size() &&
           "Invalid start position");

}

pair<double **, double **> LUDecompose(vector<vector<double> > &matrix,
                                       size_t n,
                                       optional<vector<int>> &p,
                                       optional<std::ostream> &out) {
    vector<vector<double> > L(n, vector<double>(n));
    vector<vector<double> > U(n, vector<double>(n));

    std::copy(matrix.begin(), matrix.end(), U);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = i; j < n - 1; j++) {

        }
    }

    return std::pair<double **, double **>();
}
