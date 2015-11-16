#pragma once

#include <vector>
#include <iostream>
#include <math.h>
#include <fstream>
class vect
{ 
public:
	std::vector<double> V;
	vect(void);
	vect(int n);
	void read(std::ifstream &vect);
	void make(int n);
	~vect(void);
	vect& operator=(const vect& newVect);
	double operator* (const vect& newVect);
	vect& operator-(vect& newVect);
	vect& operator+(vect& newVect); 
	double norm();
};

