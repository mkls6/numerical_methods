#ifndef LINEAR_ALGEBRA_HPP
#define LINEAR_ALGEBRA_HPP

#include <utility>
#include <vector>
#include <optional>
#include <ostream>

using std::pair, std::vector, std::optional;

pair<double **, double **> LUDecompose(vector<vector<double> > &,
                                       size_t,
                                       optional<vector<int>> &,
                                       optional<std::ostream> &);

#endif
