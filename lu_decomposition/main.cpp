#include <iostream>
#include <fstream>
#include "hpp/matrix.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);
    Matrix *m = Matrix::ReadMatrix(inputFile);

    auto [L, U, permutationsCount] = m->LUDecompose(std::nullopt);
    std::cout << "U:\n" << U << "\nL:\n" << L << "\nP_N: " << permutationsCount << std::endl;
    std::cout << Matrix::LUDeterminant(L, U, permutationsCount);

    return 0;
}
