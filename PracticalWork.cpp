#include "pch.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <fstream> 
using namespace std;
double p1(double x)
{
	return x - pow(x, 2);
}
double q1(double x)
{
	return pow(pow(sin(x), 2) + cos(x), 1 / 3);
}
double f1(double x)
{
	return pow(x, 3) - pow(x, 2) + sqrt(x);
}
double p2(double x)
{
	return sqrt(x*sin(x) + pow(x, 2));
}
double q2(double x)
{
	return sin(x) + cos(x);
}
double f2(double x)
{
	return pow(pow(pow(x, 3) - pow(x, 2), 0.5) + x, 1/3);
}
/********МЕТОД КОНЕЧНЫХ РАЗНОСТЕЙ********/
//	для уравнения вида
//	eps * y'' + p(x) * y' + q(x) * y + f(x) = 0		
//	с граничными условиями вида 
//	-alpha1 * y'(0) + alpha2 * y(0) = gamma1
//	beta1 * y'(1) + beta2 * y(1) = gamma2
void FiniteDifferenceMethod(double alpha1, double alpha2, double beta1, double beta2, double gamma1, double gamma2, double epsilon, double p(double x), double q(double x), double f(double x), const char string[])
{
	double a = 0, b = 1;
	int n = 100;
	double h = (b - a) / n;
	vector <double> x(n + 1, 0), y(n + 1, 0);
	//создаем сетку иксов
	for (int i = 0; i < n + 1; i++)
	{
		x[i] = a + i * h;
	}
	//создание векторов для значений трехдиаг. матрицы.
	vector <double> A(n, 0), B(n + 1, 0), C(n, 0), F(n + 1, 0);
	//заполняем "нижнюю" диагональ A
	for (int i = 0; i < n - 1; i++)
	{
		A[i] = epsilon - p(x[i + 1])*h / 2;
	}
	A[n - 1] = -beta1;
	//заполняем главную диагональ B
	B[0] = alpha1 + alpha2 * h;
	for (int i = 1; i < n; i++)
	{
		B[i] = pow(h, 2)* q(x[i]) - 2 * epsilon;
	}
	B[n] = beta1 + beta2 * h;
	//заполняем "верхнюю" диагональ C
	C[0] = - alpha1;
	for (int i = 1; i < n; i++)
	{
		C[i] = epsilon + p(x[i])*h / 2;
	}
	//заполняем столбец свободных членов F
	F[0] = gamma1 * h;
	for (int i = 1; i < n; i++)
	{
		F[i] = -f(x[i])* pow(h, 2);
	}
	F[n] = gamma2 * h;
	/**********МЕТОД ПРОГОНКИ**********/
	//создаем дополнительные вектора для коэффициентов P и Q
	vector <double> P(n, 0), Q(n + 1, 0);
	//заполняем первые элементы векторов P и Q
	P[0] = - C[0] / B[0];
	Q[0] = F[0] / B[0];
	//вычисляем остальные коэффициенты P
	for (int i = 1; i < n; i++)
	{
		P[i] = - C[i] / (B[i] + A[i - 1] * P[i - 1]);
	}
	//вычисляем остальные коэффициенты Q
	for (int i = 1; i < n + 1; i++)
	{
		Q[i] = (F[i] - A[i - 1] * Q[i - 1]) / (B[i] + A[i - 1] * P[i - 1]);
	}
	//вычисляем значения y
	y[n] = Q[n];
	for (int i = n - 1; i > -1; i--)
	{
		y[i] = P[i] * y[i + 1] + Q[i];
	}
	ofstream fout(string);
	for (int i = 0; i < n + 1; i++)
	{
		fout << x[i] << "," << y[i] << endl;
	}
	fout.close();
}
int main()
{
	FiniteDifferenceMethod(1, 0, 1, 0, 4, 7, 1,p1,q1,f1, "var3_1.txt");
	FiniteDifferenceMethod(1, 0, 1, 0, 4, 7, 0.1, p1, q1, f1, "var3_0.1.txt");
	FiniteDifferenceMethod(1, 0, 1, 0, 4, 7, 0.01, p1, q1, f1, "var3_0.01.txt");
	FiniteDifferenceMethod(1, 0, 1, 0, 4, 7, 0.001, p1, q1, f1, "var3_0.001.txt");

	FiniteDifferenceMethod(0, 1, 0, 1, 1, 4, 1, p2, q2, f2, "var18_1.txt");
	FiniteDifferenceMethod(0, 1, 0, 1, 1, 4, 0.1, p2, q2, f2, "var18_0.1.txt");
	FiniteDifferenceMethod(0, 1, 0, 1, 1, 4, 0.01, p2, q2, f2, "var18_0.01.txt");
	FiniteDifferenceMethod(0, 1, 0, 1, 1, 4, 0.001, p2, q2, f2, "var18_0.001.txt");

	cout << "Good work!";
}