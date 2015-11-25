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
	vect y;
	int maxiter;
	double eps;
	int iter;
	double normR;
public:
	SLAU(void);
	SLAU(std::ifstream &size, int flag); 
	SLAU( std::ifstream &size);
	~SLAU(void);
	void LUdec(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D);
	void Direct(std::vector<double> &L, std::vector<double> &D, vect &y, vect &F);
	void Reverse(std::vector<double> &U, vect &x, vect &y);
	void Direct(std::vector<double> &L, vect &y, vect &F);
	void Reverse(std::vector<double> &U, std::vector<double> &D, vect &x, vect &y);
	void LOS();
	void LOS_D();
	void LOS_LU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D);
	void MCG();
	void MCG_D();
	void MCG_LU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D);
	void output(std::ofstream &X);
	void MCG_PLUS();
	void MCG_D_PLUS();
	void LU_SLAU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D);
};

