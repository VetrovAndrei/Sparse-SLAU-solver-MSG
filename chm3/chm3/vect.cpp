#include "vect.h"


vect::vect(void)
{
}

vect::vect(int n)
{
	V.resize(n);
}

void vect::make(int n)
{
	V.resize(n);
}

vect::~vect(void)
{
}

void vect::read(std::ifstream &vect)
{
	for (int i = 0; i < V.size(); i++)
	{
		vect >> V[i];
	}
}

vect& vect::operator=(const vect& newVect)
	{
		if (this != &newVect)
			this->V = newVect.V;
		return *this;
	}

double vect::operator* (const vect& newVect)
{
	double sum = 0;
	if (this->V.size() != newVect.V.size())
		return sum;
	else
	{
		for(int i = 0; i < this->V.size(); i++)
			sum += this->V[i] * newVect.V[i];
	}
	return sum;
}

double vect::norm()
{
	double norma = 0;
	for(int i = 0; i < this->V.size(); i++)
	{
			norma += pow(this->V[i],2.0);
	}
	norma = sqrt(norma);
	return norma;
}

vect& vect::operator-(vect& newVect)
{
	if (this->V.size() != newVect.V.size())
	{
		return *this;
	}
	else
	{
		for (int i = 0; i < this->V.size(); i++)
		{
			newVect.V[i] = this->V[i] - newVect.V[i]; 
		}
		return newVect;
	}
}

vect& vect::operator+(vect& newVect)
{
	if (this->V.size() != newVect.V.size())
	{
		return *this;
	}
	else
	{
		for (int i = 0; i < this->V.size(); i++)
		{
			this->V[i] = this->V[i] + newVect.V[i]; 
		}
		return *this;
	}
}