#include <iostream>
#include <fstream>
#include <string>
#include "hpp/matrix.hpp"
#include "hpp/linear_algebra.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);
    size_t n;

    inputFile >> n;
    Matrix *m = Matrix::ReadMatrix(inputFile, n);
    auto M = *m;

    vector<double> b = {18.1000000, 60.7000000, -31.3000000, 11.9000000};
    auto B = Matrix(b);

    auto [L, U, P, permutationsCount] = M.LUDecompose();
    auto P_T = ~P;

    B = B * P_T;
    auto x = Matrix::LUSolve(L, U, B[0]);
    auto X = Matrix(x);

    std::cout << "U:\n" << U << "\nL:\n" << L << "\nP_N: " << permutationsCount << std::endl;
    std::cout << Matrix::LUDeterminant(L, U, permutationsCount) << "\n";
    std::cout << "\nComputed X:\n";
    for (auto i : x) {
        std::cout << i << " ";
    }
    std::cout << "\nL * U:\n" << L * U << "\n";

    auto inverseM = ~Matrix::LUInverseMatrix(L, U, P_T);
    std::cout << "Inverse m:\n" << inverseM << "\n";
    std::cout << "LU - PA:\n" << L * U - P * M << "\n";
    std::cout << "Ax - b:\n" << P * M * ~X - ~B << "\n";
    std::cout << "Норма 1: " << LAlgebra::CubicNorm(M) * LAlgebra::CubicNorm(inverseM) << "\n";
    std::cout << "Норма 2: " << LAlgebra::OctahedralNorm(M) * LAlgebra::OctahedralNorm(inverseM) << "\n";

    return 0;
}
