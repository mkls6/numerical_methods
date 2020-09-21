#include <iostream>
#include <fstream>
#include "hpp/matrix.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);

    Matrix *m = Matrix::ReadMatrix(inputFile);

    vector<double> b = {18.1000000, 60.7000000, -31.3000000, 11.9000000};

    auto [U, L, permutationsCount] = m->LUDecompose(&b);

    std::cout << "U:\n" << U << "\nL:\n" << L << "\nP_N: " << permutationsCount << std::endl;
    std::cout << Matrix::LUDeterminant(L, U, permutationsCount) << "\n";
    std::cout << "\nComputed X:\n";
    for (auto i : Matrix::LUSolve(L, U, b)) {
        std::cout << i << " ";
    }
    std::cout << "\n";

    return 0;
}
