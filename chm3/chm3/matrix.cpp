#include "matrix.h"


void matrix::make(std::ifstream &size)
{
	size >> n;
	std::ifstream Ggu("ggu.txt");
	std::ifstream Ggl("ggl.txt");
	std::ifstream Jg("jg.txt");
	std::ifstream Ig("ig.txt");
	std::ifstream vect("vektor.txt");
	std::ifstream Di("di.txt");
	F.make(n);
	F.read(vect);
	normF = F.norm();
	ig.resize(n + 1);
	for (int i = 0; i < n + 1; i++)
	{
		Ig >> ig[i];
	}
	if (ig[0] == 1)
	{
		for (int i = 0; i < n + 1; i++)
		{
			ig[i]--;
		}
	}
	di.resize(n);
	for (int i = 0; i < n; i++)
	{
		Di >> di[i];
	}
	ggu.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		Ggu >> ggu[i];
	}
	ggl.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		Ggl >> ggl[i];
	}
	jg.resize(ig[n]);
	for (int i = 0; i < ig[n]; i++)
	{
		Jg >> jg[i];
	}
	if (jg[0] == 1)
	{
		for (int i = 0; i < ig[n]; i++)
		{
			jg[i]--;
		}
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
	int k1,k2;
	for (int i = 0; i < n; i++)
	{
		Y.V[i] = di[i] * X.V[i];
		k1 = ig[i];
		k2 = ig[i+1];
		for(int k = k1; k < k2; k++)
		{
			int j = jg[k];
			Y.V[i] += ggl[k] * X.V[j];
			Y.V[j] += ggu[k] * X.V[i];
		}
	}
}

// умножение транспонированной матрицы на вектор. 1: вектор 2: результат
void matrix::multMTV(vect &X, vect &Y)
{
	int k1, k2;
	for (int i = 0; i < n; i++)
	{
		Y.V[i] = di[i] * X.V[i];
		k1 = ig[i];
		k2 = ig[i+1];
		for(int k = k1; k < k2; k++)
		{
			int j = jg[k];
			Y.V[i] += ggu[k] * X.V[j];
			Y.V[j] += ggl[k] * X.V[i];
		}
	}
}

//деление вектора на диагональ
void matrix::divDi(vect &X)
{
	for (int i = 0; i < n; i++)
	{
		X.V[i] /= di[i];
	}
}

void matrix::Hilbert(std::ifstream &size)
{
	size >> n;
	ig.resize(n + 1);
	ig[0] = 0;
	for (int i = 1; i < n + 1; i++)
	{
		ig[i] = ig[i-1] + i - 1;
	}
	di.resize(n);
	ggu.resize(ig[n]);
	ggl.resize(ig[n]);
	jg.resize(ig[n]);
	int k1, k2;
	for (int i = 0; i < n; i++)
	{
		di[i] = double(1.0 / (2.0 * i + 1.0));
		k1 = ig[i];
		k2 = ig[i+1];
		int l = 0;
		for(int k = k1; k < k2; k++, l++)
		{
			ggl[k] = double(1.0 / (i + l + 1.0));
			ggu[k] = ggl[k];
			jg[k] = l;
		}
	}
	vect x(n);
	for (int i = 0; i < n; i++)
	{
		x.V[i] = i + 1;
	}
	F.make(n);
	this->multMV(x,F);
	normF = F.norm();
}