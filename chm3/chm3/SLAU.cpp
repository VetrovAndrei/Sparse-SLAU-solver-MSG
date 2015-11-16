#include "SLAU.h"


SLAU::SLAU(void)
{
}

SLAU::SLAU(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size, std::ifstream& X)
{
	A.make(matrix, vect, size);
	size >> maxiter >> eps;
	x0.make(A.n);
	x0.read(X);
	r.make(A.n);
	z.make(A.n);
	p.make(A.n);
	xK.make(A.n);
	Ar.make(A.n);
}

SLAU::~SLAU(void)
{
}

void SLAU::LOS()
{
	double scalPP = 0,scalPR = 0, scalRR = 0, scalPAr = 0;
	double a = 0, b = 0;
	int iter;
	bool exit = 1;
	A.multMV(x0,r);
	A.F - r;
	z = r;
	A.multMV(z,p);
	scalRR = r * r;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( abs(scalRR) < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		scalPP = p * p;
		scalPR = p * r;
			a = scalPR/scalPP;
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - p.V[i] * a;
		}
		A.multMV(r, Ar);
		scalPAr = p * Ar;
		b = -(scalPAr/scalPP);
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
			p.V[i] = Ar.V[i] + p.V[i] * b;
		}	
		scalRR = scalRR - pow(a,2.0) * scalPP;
	}
}