#include "Random.h"
#include "Matrix.h"
#include "Timer.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>

#include "qr_algorithm.h"
#include "power_method.h"
#include "nonlinear_equation.h"


#define PI           3.14159265358979323846  
#define E			 2.71828182845904523536



void time_pm()
{
	std::ofstream of("output\\time_power_method.txt");
	Timer t;
	for (int i = 3; i <= 2000; i += 10)
	{
		auto k = complex::Matrix::Generate(i, i, default_complex_generator(-1, 1));

		t.start();
		power_method::base(k);
		t.stop();

		of << t.elapsed() << "\n";
	}

	of.close();
}

void time_qr()
{
	std::ofstream of("output\\time_qr_algorithm.txt");
	Timer t;
	for (int i = 3; i <= 2000; i += 10)
	{
		auto k = Matrix::Generate(i, i, default_double_generator(-1, 1));

		t.start();
		QR(k);
		t.stop();

		of << t.elapsed() << "\n";
	}

	of.close();
}


void pm()
{
	std::ofstream of1("output\\pm_a1.txt"),
				  of2("output\\pm_a2.txt");
	power_method::base(complex::Matrix::FromFile("input\\A1.txt"), of1, true);
	power_method::base(complex::Matrix::FromFile("input\\A2.txt"), of2, true);

	of1.close();
	of2.close();
}

void qr()
{
	std::ofstream of1("output\\qr_a1.txt"),
				  of2("output\\qr_a2.txt");
	QR(Matrix::FromFile("input\\A1.txt"), of1);
	QR(Matrix::FromFile("input\\A2.txt"), of2);

	of1.close();
	of2.close();
}

double f(const double& x)
{
	double y = std::pow(x, 9) + PI;
	y *= std::cos(std::log(x * x + 1));
	y /= std::pow(E, x * x);
	y -= x / 2022;
	return y;
}

double df(const double& x)
{
	double y = 2 * std::pow(E, -(x * x)) * x * (std::pow(x, 9) + PI) * std::sin(std::log(x * x + 1));
	y = -y;
	y /= x * x + 1;
	double t = std::pow(E, -(x * x)) * x * (2 * std::pow(x, 9) - 9 * std::pow(x, 7) + 2 * PI) * std::cos(std::log(x * x + 1)) - 1 / 2022;
	y -= t;
	return y;
}

void ne(double&& from, double&& to)
{
	std::cout << "\nINTERVAL: [" << from << "; " << to << "]";
	
	int i = 0;

	i = bisexion(f, from, to, 1e-4);

	std::cout << "\nBISEXION: [" << from << "; " << to << "]";
	std::cout << "\n" << i << " ITERATIONS";

	double x = from;

	auto i1 = newton(f, df, from, 1e-15);
	auto i2 = discrete_newton(f, x, 0.0001, 1e-15);

	
	std::cout << "\n\nDerivative res: f(" << from << ") =" << f(from);
	std::cout << "\n" << i1 << " ITERATIONS";
	std::cout << "\nDiscrete   res: f(" << x    << ") =" << f(x);
	std::cout << "\n" << i2 << " ITERATIONS";
	std::cout << "\n\n-----------------------------------------------------------\n";
}

void linear()
{
	ne(-2.5, -1.7);
	ne(-1.4, -0.5);
	ne(1.57, 2.4);
}

int main()
{
	//pm();
	//qr();

	linear();


	return 0;
}