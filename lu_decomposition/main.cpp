#include <iostream>
#include <fstream>
#include "hpp/matrix.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);

    Matrix *m = Matrix::ReadMatrix(inputFile);
    auto M = *m;

    vector<double> b = {18.1000000, 60.7000000, -31.3000000, 11.9000000};
    vector<vector<double> > t;
    t.push_back(b);
    auto B = Matrix(t);

    auto [L, U, P, permutationsCount] = M.LUDecompose();
    auto P_T = ~P;

    B = B * P_T;
    auto x = Matrix::LUSolve(L, U, B[0]);
    vector<vector<double>> tx; tx.push_back(x);
    auto X = Matrix(tx);

    std::cout << "U:\n" << U << "\nL:\n" << L << "\nP_N: " << permutationsCount << std::endl;
    std::cout << Matrix::LUDeterminant(L, U, permutationsCount) << "\n";
    std::cout << "\nComputed X:\n";
    for (auto i : x) {
        std::cout << i << " ";
    }
    std::cout << "\nL * U:\n" << L * U << "\n";
    std::cout << "Inverse m:\n" << ~Matrix::LUInverseMatrix(L, U, P_T) << "\n";
    std::cout << "LU - PA:\n" << L * U - P * M << "\n";
    std::cout << "Ax - b:\n" << P * M * ~X - ~B << "\n";

    return 0;
}
