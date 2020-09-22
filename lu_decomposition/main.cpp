#include <iostream>
#include <fstream>
#include <string>
#include "hpp/matrix.hpp"
#include "hpp/linear_algebra.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    if (filePath.empty()) {
        std::cout << "Will try to open input.txt\n";
        filePath = "input.txt";
    }

    std::ifstream inputFile;

    try {
        inputFile.open(filePath);
    } catch (std::exception &e) {
        std::cout << "Failed to open file " << filePath << "." << std::endl;
        exit(1);
    }
    size_t n;

    inputFile >> n;
    Matrix *a = Matrix::ReadMatrix(inputFile, n);
    auto A = *a;

    vector<double> b(n, 0);
    for (size_t i = 0; i < n; i++)
        inputFile >> b[i];

    vector<double> referenceX(n, 0);
    for (size_t i = 0; i < n; i++)
        inputFile >> referenceX[i];

    auto B = Matrix(b);

    auto [L, U, P, permutationsCount] = A.LUDecompose();
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
    std::cout << "\nX - X_ref:\n";
    for (size_t i = 0; i < n; i++) {
        std::cout << x[i] - referenceX[i] << " ";
    }
    std::cout << "\nL * U:\n" << L * U << "\n";

    auto inverseM = Matrix::LUInverseMatrix(L, U, P_T);
    std::cout << "Inverse a:\n" << inverseM << "\n";
    std::cout << "LU - PA:\n" << L * U - P * A << "\n";
    std::cout << "Ax - b:\n" << P * A * ~X - ~B << "\n";
    std::cout << "Cubic norm: " << LAlgebra::CubicNorm(A) * LAlgebra::CubicNorm(inverseM) << "\n";
    std::cout << "Octahedral norm: " << LAlgebra::OctahedralNorm(A) * LAlgebra::OctahedralNorm(inverseM) << "\n";
    std::cout << "Euclidean norm: " << LAlgebra::EuclideanNorm(A) * LAlgebra::EuclideanNorm(inverseM) << "\n";

    return 0;
}
