#include "matrix.h"


void matrix::make(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size)
{
	size >> n;
	F.make(n);
	F.read(vect);
	normF = F.norm();
	ig.resize(n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		matrix >> ig[i];
	}
	di.resize(n);
	for (int i = 0; i < n; i++)
	{
		matrix >> di[i];
	}
	ggu.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		matrix >> ggu[i];
	}
	ggl.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		matrix >> ggl[i];
	}
	jg.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		matrix >> jg[i];
	}
}

matrix::matrix(void)
{
}

matrix::~matrix(void)
{
}

// умножение матрицы на вектор. 1: вектор 2: результат
void matrix::multMV(vect &X, vect &Y)
{
	// инициализация результата
	for (int i = 0; i < n; i++)
	{
		Y.V[i] = di[i] * X.V[i];
	}
	// проход по внедиагональным элементам
	for (int i = 0; i < n; i++)
	{
		for(int k = ig[i]; k < ig[i+1]; k++)
		{
			Y.V[i] += ggl[k] * X.V[jg[k]];
			Y.V[jg[k]] += ggu[k] * X.V[i];
		}
	}
}