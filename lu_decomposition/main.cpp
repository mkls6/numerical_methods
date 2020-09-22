#include <iostream>
#include <fstream>
#include "hpp/matrix.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);

    Matrix *m = Matrix::ReadMatrix(inputFile);

    vector<double> b = {18.1000000, 60.7000000, -31.3000000, 11.9000000};
    vector<vector<double> > t;
    t.push_back(b);
    auto B = Matrix(t);

    auto [U, L, P, permutationsCount] = m->LUDecompose();

    P = ~P;
    B = B * P;
    *m = *m * P;

    std::cout << "U:\n" << U << "\nL:\n" << L << "\nP_N: " << permutationsCount << std::endl;
    std::cout << Matrix::LUDeterminant(L, U, permutationsCount) << "\n";
    std::cout << "\nComputed X:\n";
    for (auto i : Matrix::LUSolve(L, U, B[0])) {
        std::cout << i << " ";
    }
    std::cout << "\nL * U:\n" << L * U << "\n";
    std::cout << "Inverse m:\n" << ~Matrix::LUInverseMatrix(L, U, P) << "\n";

    return 0;
}
