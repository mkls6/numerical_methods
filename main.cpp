#include <iostream>
#include <fstream>
#include "hpp/matrix.hpp"

int main() {
    std::string filePath;
    std::cin >> filePath;
    std::ifstream inputFile(filePath);
    Matrix *m = Matrix::ReadMatrix(inputFile);

    auto res = m->LUDecompose(std::nullopt);
    std::cout << "U:\n" << res.first << "\nL:\n" << res.second << std::endl;

    return 0;
}
