#include "SLAU.h"

void main()
{
	std::ifstream size ("size.txt");
	std::ofstream X("X.txt");
	SLAU mat(size, 1);
	std::vector<double> L;
	std::vector<double> U;
	std::vector<double> D;
	mat.LUdec(L, U, D);
	//mat.LU_SLAU(L, U, D);
	mat.MCG_LU(L, U, D);
	mat.output(X);
	system("pause");
}