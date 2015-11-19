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
	Ar.make(A.n);
}

SLAU::~SLAU(void)
{
}

//ЛОС без предобусловливания
void SLAU::LOS()
{
	double scalPP = 0,scalPR = 0, scalRR = 0, scalPAr = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r);
	A.F - r;
	z = r;
	A.multMV(z,p);
	scalRR = r * r;
	normR = sqrt(scalRR) / A.normF; 
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
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
		normR = sqrt(scalRR)/A.normF;
	}
}

//ЛОС с предоубсловливанием диагонали
void SLAU::LOS_D()
{
	double scalPP = 0,scalPR = 0, scalRR = 0, scalPAr = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r);
	A.F - r;
	A.divDi(r);
	z = r;
	A.multMV(z,p);
	A.divDi(p);
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
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
		A.divDi(Ar);
		scalPAr = p * Ar;
		b = -(scalPAr/scalPP);
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
			p.V[i] = Ar.V[i] + p.V[i] * b;
		}	
		scalRR = scalRR - pow(a,2.0) * scalPP;
		normR = sqrt(scalRR)/A.normF;
	}
}


void SLAU::MCG()
{
	double scalRR = 0, scalAzZ = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r); // R = A*x
	A.F - r; // R = F - R
	A.multMTV(r,z); // z = At*R
	r = z; // r = z
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		A.multMV(z,p); // P = A*Z
		A.multMTV(p,Ar); // Ar = At*P
		scalAzZ = Ar * z; //(Ar,z)
		a = scalRR/scalAzZ; 
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a;
		}
		b = - 1.0/scalRR;
		scalRR = r * r;
		b *= scalRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
		}	
		normR = sqrt(scalRR)/A.normF;
	}
}

void SLAU::MCG_D()
{
	double scalRR = 0, scalAzZ = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r); //r = A*x
	A.F - r; // r = F-r
	A.multMTV(r,z);// z = At*r
	A.divDi(z); // z = A(-1)*z
	A.divDi(z); // z = A(-1)*z
	r = z; // r = z;
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		A.multMV(z,p); // p = A*z
		A.multMTV(p,Ar); // Ar = At*p
		A.divDi(Ar);
		A.divDi(Ar);
		scalAzZ = Ar * z; // (Ar,z)
		a = scalRR/scalAzZ;
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a; // Ar = D-2 * At*A*z
		}
		b = - 1.0/scalRR;
		scalRR = r * r;
		b *= scalRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = p.V[i] + z.V[i] * b;
		}	
		normR = sqrt(scalRR)/A.normF;
	}
}

void SLAU::LUdec(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D)
{

}


void SLAU::output(std::ofstream &X)
{
	X << iter << std::endl;
	X << normR << std::endl;
	for (int i = 0; i < A.n; i++)
	{
		X << x0.V[i] << std::endl;
	}
}