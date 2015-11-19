#pragma once
#include "vect.h"
class matrix
{
public:
	std::vector<double> di; //диагональные элементы
	std::vector<double> ggl; // элементы нижней диагонали
	std::vector<double> ggu; // элементы верхней диагонали
	vect F; // вектор правой части
	std::vector<int> jg; // массив с индексами начала строк(столбцов) в массиве ggl(ggu) 
	std::vector<int> ig; // номера столбцов(строк) элементов в массиве ggl(ggu)
	int n; // размерность матрицы
	double normF; // норма вектора F
	matrix(void);
	void make(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size);
	~matrix(void);
	void multMV(vect &X, vect &Y);
	void multMTV(vect &X, vect &Y);
	void divDi(vect &X);
	void Hilbert(std::ifstream &size);
};

