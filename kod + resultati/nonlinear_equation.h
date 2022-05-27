#pragma once

#include <functional>

double bisexion(std::function<double(const double&)> f, double& from, double& to, const double& error)
{
	bool negative_first = f(from) < 0;

	int i = 0;
	while (to - from > error)
	{
		double x = (from + to) / 2;
		double y = f(x);
		if (y <= 0)
		{
			if (negative_first)
				from = x;
			else
				to = x;
		}
		else
		{
			if (negative_first)
				to = x;
			else
				from = x;
		}

		++i;
	}

	return i;
}

double newton(std::function<double(const double&)> function, decltype(function) derivative, double& x, const double& error)
{
	double x2 = x + 5;

	int i = 0;
	while (std::abs(x2 - x) > error)
	{
		double y = function(x);
		double d = derivative(x);
		x2 = x;
		x = x2 - y / d;

		++i;
	}

	return i;
}

double discrete_newton(std::function<double(const double&)> function, double& x, const double& delta, const double& error)
{
	double x2 = x + 5;

	int i = 0;
	while (std::abs(x2 - x) > error)
	{
		double y = function(x);
		double d = function(x + delta) / delta;
		x2 = x;
		x = x2 - y / d;

		++i;
	}

	return i;
}