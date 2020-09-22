#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include "hpp/matrix.hpp"
#include "hpp/linear_algebra.hpp"

int main() {
    std::string filePath = "input.txt";
    std::ifstream inputFile;
    std::ofstream outputFile;
    size_t n;
    Matrix *a;
    vector<double> b;
    vector<double> referenceX;

    inputFile.exceptions(std::ifstream::badbit | std::ifstream::failbit);
    outputFile.exceptions(std::ofstream::badbit | std::ofstream::failbit);

    try {
        inputFile.open(filePath);
    } catch (std::ios_base::failure &e) {
        std::cout << "Failed to open file " << filePath << ".\n"
                  << "Exception message:\n" << e.what() << std::endl;
        exit(1);
    }

    try {
        inputFile >> n;
        a = Matrix::ReadMatrix(inputFile, n);
    } catch (std::ios_base::failure &) {
        std::cout << "Failed to read matrix from input file." << std::endl;
        exit(1);
    }
    auto A = *a;

    b.resize(n, 0);
    try {
        for (size_t i = 0; i < n; i++)
            inputFile >> b[i];
    } catch (std::ios_base::failure &) {
        std::cout << "Failed to read b vector from file." << std::endl;
        exit(1);
    }

    referenceX.resize(n, 0);
    try {
        for (size_t i = 0; i < n; i++)
            inputFile >> referenceX[i];
    } catch (std::ios_base::failure &) {
        std::cout << "Failed to read X vector from input file." << std::endl;
        exit(1);
    }
    inputFile.close();

    try {
        outputFile.open("output.txt");
        std::cout.rdbuf(outputFile.rdbuf());
    } catch (std::ios_base::failure &) {
        std::cout << "Failed to write to output file. Writing to stdout instead.\n";
    }

    auto B = Matrix(b);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Matrix A:\n" << A << "\n";

    auto [L, U, P, permutationsCount] = A.LUDecompose();
    auto P_T = ~P;

    B = B * P_T;
    auto x = Matrix::LUSolve(L, U, B[0]);
    auto X = Matrix(x);
    auto inverseA = Matrix::LUInverseMatrix(L, U, P_T);

    std::cout << "Determinant = " << Matrix::LUDeterminant(L, U, permutationsCount) << "\n";

    std::cout << "\nComputed X:\n";
    for (auto i : x) {
        std::cout << std::setw(16) << i;
    }

    std::cout << "\n\nX - X_ref:\n" << std::scientific << std::setprecision(6);
    for (size_t i = 0; i < n; i++) {
        std::cout << std::setw(16) << x[i] - referenceX[i];
    }

    std::cout << "\n\nP * L * U:\n" << std::fixed << std::setprecision(8) << P_T * L * U << "\n";

    std::cout << "Inverse A:\n" << inverseA << "\n";
    std::cout << "LU - PA:\n" << std::scientific << std::setprecision(6) << L * U - P * A << "\n";
    std::cout << "Ax - b:\n" << P * A * ~X - ~B << "\n";
    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Cubic norm: " << LAlgebra::CubicNorm(A) * LAlgebra::CubicNorm(inverseA) << "\n";
    std::cout << "Octahedral norm: " << LAlgebra::OctahedralNorm(A) * LAlgebra::OctahedralNorm(inverseA) << "\n";
    std::cout << "Euclidean norm: " << LAlgebra::EuclideanNorm(A) * LAlgebra::EuclideanNorm(inverseA) << "\n";

    if (outputFile.is_open())
        outputFile.close();

    return 0;
}
