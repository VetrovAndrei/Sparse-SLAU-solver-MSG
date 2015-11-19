#pragma once
#include "matrix.h"
class SLAU
{
private:
	matrix A;
	vect p;
	vect z;
	vect r;
	vect x0;
	vect Ar;
	int maxiter;
	double eps;
	int iter;
	double normR;
public:
	SLAU(void);
	SLAU(std::ifstream &size, std::ifstream &X); 
	SLAU(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size, std::ifstream& X);
	~SLAU(void);
	void LUdec(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D);
	void LOS();
	void LOS_D();
	void MCG();
	void MCG_D();
	void output(std::ofstream &X);
};

