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
	vect xK;
	vect Ar;
	int maxiter;
	double eps;

public:
	SLAU(void);
	SLAU(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size, std::ifstream& X);
	~SLAU(void);
	void LOS();
};

