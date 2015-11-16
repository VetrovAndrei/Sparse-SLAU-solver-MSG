#include "SLAU.h"

void main()
{
	std::ifstream size ("size.txt");
	std::ifstream matr("matrix.txt");
	std::ifstream vect("vector.txt");
	std::ifstream nachX("x0.txt");
	std::ofstream X("X.txt");
	SLAU mat(matr, vect, size, nachX);
	mat.LOS();
	system("pause");
}