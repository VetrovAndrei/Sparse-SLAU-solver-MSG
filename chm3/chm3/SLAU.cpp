#include "SLAU.h"


SLAU::SLAU(void)
{
}

SLAU::SLAU(std::ifstream &size, int flag)
{
	A.make(size);
	size >> maxiter >> eps;
	x0.make(A.n);
	r.make(A.n);
	z.make(A.n);
	p.make(A.n);
	Ar.make(A.n);
	y.make(A.n);
}

SLAU::SLAU(std::ifstream &size)
{
	A.Hilbert(size);
	size >> maxiter >> eps;
	x0.make(A.n);
	r.make(A.n);
	z.make(A.n);
	p.make(A.n);
	Ar.make(A.n);
	y.make(A.n);
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
		//scalRR = scalRR - pow(a,2.0) * scalPP;
		scalRR = r * r;
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
		//scalRR = scalRR - pow(a,2.0) * scalPP;
		scalRR = r * r;
		normR = sqrt(scalRR)/A.normF;
	}
}

void SLAU::LOS_LU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D)
{
	double scalPP = 0,scalPR = 0, scalRR = 0, scalPAr = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,y); // y = A * x
	A.F - y;		// y = F - y
	Direct(L,D,r,y);// r = L-1 * y
	Reverse(U,z,r); // z = U-1 * r
	A.multMV(z,y);  // y = A * z
	Direct(L,D,p,y);// p = L-1 * y 
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
		Reverse(U, y, r);
		A.multMV(y, Ar);
		Direct(L, D, Ar, Ar);
		scalPAr = p * Ar;
		b = -(scalPAr/scalPP);
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = y.V[i] + z.V[i] * b;
			p.V[i] = Ar.V[i] + p.V[i] * b;
		}	
		//scalRR = scalRR - pow(a,2.0) * scalPP;
		scalRR = r * r;
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
		b =  1.0/scalRR;
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
	A.divDi(r); // r = A(-1)*r
	A.divDi(r); // r = A(-1)*r
	A.multMTV(r,z);// z = At*r
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
		A.divDi(p);
		A.divDi(p);
		A.multMTV(p,Ar); // Ar = At*p
		scalAzZ = Ar * z; // (Ar,z)
		a = scalRR/scalAzZ;
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a; // Ar = D-2 * At*A*z
		}
		b = 1.0/scalRR;
		scalRR = r * r;
		b *= scalRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
		}	
		normR = sqrt(scalRR)/A.normF;
	}
}

void SLAU::MCG_LU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D)
{
	double scalRR = 0, scalAzZ = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r); // r = A*x
	A.F - r; // R = F - R
	Direct(L, D, r, r); // r = L-1 * r
	Reverse(L, D, r, r); // r = L-t * r
	A.multMTV(r,y); // y = At * r
	Direct(U, r, y); // r = U-t * y
	z = r; // r = z
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		Reverse(U, y, z); // y = U-1 * z
		A.multMV(y,p); // P = A * y
		Direct(L, D, p, p); // p = L-1 * p
		Reverse(L, D, p, p); // p = L-t * p
		A.multMTV(p,Ar); // Ar = At * P
		Direct(U, Ar, Ar); // Ar = U-t * Ar
		scalAzZ = Ar * z; //(Ar,z)
		a = scalRR/scalAzZ; 
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a;
		}
		b =  1.0/scalRR;
		scalRR = r * r;
		b *= scalRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
		}	
		normR = sqrt(scalRR)/A.normF;
	}
	Reverse(U, x0, x0); // x0 = U-1 * x0
}

void SLAU::MCG_D_PLUS()
{
	double scalRR = 0, scalAzZ = 0, scalMRR = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r); // R = A*x
	A.F - r; // R = F - R
	z = r;
	A.divDi(z);
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		A.multMV(z,Ar); // P = A*Z
		//A.multMTV(p,Ar); // Ar = At*P
		scalAzZ = Ar * z; //(Ar,z)
		p = r;
		A.divDi(p);
		scalMRR = p * r;
		a = scalMRR/scalAzZ; 
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a;
		}
		b =  1.0/scalMRR;
		p = r;
		A.divDi(p);
		scalMRR = p * r;
		b *= scalMRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = p.V[i] + z.V[i] * b;
		}	

		scalRR = r * r;
		normR = sqrt(scalRR)/A.normF;
	}
}

void SLAU::MCG_PLUS()
{
	double scalRR = 0, scalAzZ = 0;
	double a = 0, b = 0;
	bool exit = 1;
	A.multMV(x0,r); // R = A*x
	A.F - r; // R = F - R
	z = r;
	scalRR = r * r;
	normR = sqrt(scalRR)/A.normF;
	for (iter = 1; iter < maxiter && exit != 0; iter++)
	{
		if( normR < eps || iter == maxiter - 1)
		{
			exit = 0;
			break;
		}
		A.multMV(z,Ar); // P = A*Z
		//A.multMTV(p,Ar); // Ar = At*P
		scalAzZ = Ar * z; //(Ar,z)
		a = scalRR/scalAzZ; 
		for (int i = 0; i < A.n; i++)
		{
			x0.V[i] = x0.V[i] + z.V[i] * a;
			r.V[i] = r.V[i] - Ar.V[i] * a;
		}
		b =  1.0/scalRR;
		scalRR = r * r;
		b *= scalRR;
		for (int i = 0; i < A.n; i++)
		{
			z.V[i] = r.V[i] + z.V[i] * b;
		}	
		normR = sqrt(scalRR)/A.normF;
	}
}


void SLAU::LUdec(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D)
{
	L = A.ggl;
	U = A.ggu;
	D = A.di;
	double l, u, d;
	for(int k = 0; k < A.n; k++)
	{
		d = 0;
		int i0 = A.ig[k], i1 = A.ig[k + 1];
		int i = i0;
		for(; i0 < i1; i0++)
		{
			l = 0;
			u = 0;
			int j0 = i, j1 = i0;
			for(; j0 < j1; j0++)
			{
				int t0 = A.ig[A.jg[i0]], t1 = A.ig[A.jg[i0]+1];
				for (;   t0 < t1; t0++)
				{
					if (A.jg[j0] == A.jg[t0])
					{
						l += L[j0] * U[t0];
						u += L[t0] * U[j0];
					}
				}
			}
			L[i0] -= l;
			U[i0] -= u;
			U[i0] /= D[A.jg[i0]];
			d += L[i0] * U[i0];
		}
		D[k] -= d;
	}
}

// прямой ход Ly=F
void SLAU::Direct(std::vector<double> &L,std::vector<double> &D, vect &y, vect &F)
{
	y = F;
	for (int i = 0; i < A.n; i++)
	{
		double sum = 0;
		int k0 = A.ig[i], k1 = A.ig[i + 1];
		int j;
		for (; k0 < k1; k0++)
		{
			j = A.jg[k0];
			sum += y.V[j] * L[k0];
		}
		double buf = y.V[i] - sum;
		y.V[i] = buf / D[i];
	}
}

// обратный ход Ux=y
void SLAU::Reverse(std::vector<double> &U, vect &x, vect &y) 
{
	x = y;
	for(int i = A.n - 1; i >= 0; i--)
	{
		int k0 = A.ig[i], k1 = A.ig[i + 1];
		int j;
		for (; k0 < k1; k0++)
		{
			j = A.jg[k0];
			x.V[j] -= x.V[i]*U[k0];
		}
	}
}

 // прямой ход Ly=F, если на диагонали единицы
void SLAU::Direct(std::vector<double> &L, vect &y, vect &F)
{
	y = F;
	for (int i = 0; i < A.n; i++)
	{
		double sum = 0;
		int k0 = A.ig[i], k1 = A.ig[i + 1];
		int j;
		for (; k0 < k1; k0++)
		{
			j = A.jg[k0];
			sum += y.V[j] * L[k0];
		}
		y.V[i] -= sum;
	}
}
// обратный ход Ux=y, если диагональ хранится в U
void SLAU::Reverse(std::vector<double> &U, std::vector<double> &D, vect &x, vect &y) 
{
	x = y;
	for(int i = A.n - 1; i >= 0; i--)
	{
		int k0 = A.ig[i], k1 = A.ig[i + 1];
		int j;
		x.V[i] /= D[i];
		for (; k0 < k1; k0++)
		{
			j = A.jg[k0];
			x.V[j] -= x.V[i]*U[k0];
		}
		
	}
}


void SLAU::output(std::ofstream &X)
{
	X << iter << std::endl;
	X.precision(17);
	X << normR << std::endl;
	for (int i = 0; i < A.n; i++)
	{
		X << x0.V[i] << std::endl;
	}
}

void SLAU::LU_SLAU(std::vector<double> &L, std::vector<double> &U, std::vector<double> &D)
{
	y = A.F;
	Direct(L, D, y, y);
	Reverse(U, x0, y);
}