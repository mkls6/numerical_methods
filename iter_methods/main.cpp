#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include "linear_algebra.hpp""
#include "matrix.hpp"

const double eps = 1e-4;

int main()
{
	std::ifstream in("input.txt");
	std::ofstream out("out.txt");

	auto A = Matrix::ReadMatrix(in, 4);
	auto b = Matrix::ReadMatrix(in, 1, 4);
	out << "Matrix A:\n" << *A << "\nVector b:\n" << *b;
	out << "\nEuclideanNorm(A): " << LAlgebra::EuclideanNorm(*A);
	
	out << "\nSimple iteration method:\n|"<< std::setw(5) << "iter" << "|" << std::setw(14) << "Tau" << "|" << std::setw(14) << "q" << "|" <<
		std::setw(14) << "residual norm" << "|" << std::setw(14) << "criterion" << "|" << std::setw(34) << "current x\n";
	unsigned int simple_time = clock();
	auto answ = LAlgebra::SimpleIterationMethod(*A, *b, eps, out);
	simple_time = clock() - simple_time;
	
	out << "\nFastest Gradient Descent Method:\n|" << std::setw(5) << "iter" << "|" << std::setw(14) << "Tau" << "|" << std::setw(14) << "q" << "|" <<
		std::setw(14) << "residual norm" << "|" << std::setw(14) << "criterion" << "|" << std::setw(34) << "current x\n";
	unsigned int fast_grad_time = clock();
	auto answ_grad = LAlgebra::FastestGradientDescentMethod(*A, *b, eps, out);
	fast_grad_time = clock() - fast_grad_time;
	
	out << "\nSOR method:\n";
	unsigned int sor_find_time = clock();
	auto[optim_w, min_iter] = LAlgebra::SorFind(*A, *b, 1e-2, out);
	sor_find_time = clock() - sor_find_time;
	out << "\nw*= " << optim_w << " min_iter= " << min_iter << "\n";
	
	out << "\nSOR Method:\n|" << std::setw(5) << "iter" << "|" << std::setw(14) << "W" << "|" << std::setw(14) << "q" << "|" <<
		std::setw(14) << "residual norm" << "|" << std::setw(14) << "criterion" << "|" << std::setw(34) << "current x\n";
	unsigned int sor_time = clock();
	auto answ_sor = LAlgebra::Sor(*A, *b, eps, optim_w, out);
	sor_time = clock() - sor_time;

	out << "\nConjugate Gradient Method:\n|" << std::setw(5) << "iter" << "|" << std::setw(14) << "Tau" << "|" << std::setw(14) << "alpha" << "| "<< 
		std::setw(13) << "q" << "|" << std::setw(14) << "residual norm" << "|" << std::setw(14) << "criterion" << "|" << std::setw(34) << "current x\n";
	unsigned int conj_time = clock();
	auto answ_conj = LAlgebra::ConjugateGradientMethod(*A, *b, eps, out);
	conj_time = clock() - conj_time;

	out << "\nTime (ms):\n" << "Simple Iteration Method: " << simple_time << "\n" <<
			"Fastest Gradient Descent Method: " << fast_grad_time << "\n" << 
			"SOR Find Params Time: " << sor_find_time << "\n" <<
			"SOR method: " << sor_time << "\n" << 
			"Conjugate Gradient Method: " << conj_time << "\n";
	
	auto[L, U, P, s] = A->LUDecompose();
	auto inverseA = A->LUInverseMatrix(L, U, P);
	auto condA = LAlgebra::EuclideanNorm(*A) * LAlgebra::EuclideanNorm(inverseA);
	
	out << "\nCond(A) = " << condA << "\n";
	out << "\nTheoretical iteration count:\n" <<
		"Simple Iteration Method and Fastest Gradient Descent Method: " << (int)(log(1 / eps) / 2 * condA) << "\n" <<
		"SOR method: " << (int)(log(1 / eps) / 4 * sqrt(condA)) << "\n" << 
		"Conjugate Gradient Method: " << (int)(log(2 / eps) / 2 * sqrt(condA)) << "\n";

	
	auto b_vec = vector<double>(A->Columns());
	for (size_t i = 0; i < A->Columns(); i++)
		b_vec[i] = (*b)[0][i];
	auto answ_lu = A->LUSolve(L, U, b_vec);
	out << "\nAnswers:\n" <<
		std::setw(40) << "Simple Iteration Method: " << answ <<
		std::setw(40) << "Fastest Gradient Descent Method: " << answ_grad <<
		std::setw(40) << "SOR method: " << answ_sor <<
		std::setw(40) << "Conjugate Gradient Method: " << answ_conj <<
		std::setw(40) << "LU Decomposition :";
	for (size_t i = 0; i < A->Columns(); i++)
		out << std::setw(10) << answ_lu[i] << " ";

	return 0;
}

