#pragma once
#include "vect.h"
class matrix
{
public:
	std::vector<double> di; //������������ ��������
	std::vector<double> ggl; // �������� ������ ���������
	std::vector<double> ggu; // �������� ������� ���������
	vect F; // ������ ������ �����
	std::vector<int> jg; // ������ � ��������� ������ �����(��������) � ������� ggl(ggu) 
	std::vector<int> ig; // ������ ��������(�����) ��������� � ������� ggl(ggu)
	int n; // ����������� �������
	double normF; // ����� ������� F
	matrix(void);
	void make(std::ifstream &matrix,  std::ifstream &vect, std::ifstream &size);
	~matrix(void);
	void multMV(vect &X, vect &Y);
	void multMTV(vect &X, vect &Y);
	void divDi(vect &X);
	void Hilbert(std::ifstream &size);
};

