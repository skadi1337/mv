#include "Task.h"

const double				  Task::MAX_RANDOM			= std::pow(2, 2.5);
const double				  Task::MIN_RANDOM			= -Task::MAX_RANDOM;

const std::function<double()> Task::RANDOM_GENERATOR	= default_double_generator(Task::MIN_RANDOM, Task::MAX_RANDOM);


Task::Task(size_t rows, size_t columns, std::function<double()> generator)
	: _rows		(rows)
	, _columns	(columns)
	, _generator (generator)
{ }

void Task::task1()
{
	auto data = std::vector<std::vector<double>>(_rows, std::vector<double>(_columns));

	double sum = 0, 
		   abs;
	std::vector<double> sums(_rows, 0);

	for (size_t i = 0; i < _rows; ++i)
	{
		sum = 0;
		for (size_t j = i + 1; j < _columns; ++j)
		{
			data[i][j] = _generator();
			data[j][i] = data[i][j];

			abs		= std::abs(data[i][j]);
			sum		+= abs;
			sums[j] += abs;
		}

		data[i][i] = 1 + sum + sums[i];
	}

	_matrix = Matrix(_rows, _columns, data);
}

void Task::task2()
{
	auto data = std::vector<std::vector<double>>(_rows, std::vector<double>(1));

	for (size_t i = 0; i < _rows; ++i)
	{
		data[i][0] = _generator();
	}

	_x = Matrix(_rows, 1, data);
	_b = _matrix.Multiply(_x);
}

void Task::task3_1()
{
	_inverse = _matrix.Inverse();
}

void Task::task3_2()
{
	_condition_number = _inverse.norm() * _matrix.norm();
}

void Task::task4()
{
	_gauss_solution = _matrix.SLAE(_b);
}

void Task::task5_1()
{
	_lup = _matrix.LUP();
}

void Task::task5_2()
{
	_lup_solution = Matrix::SLAE(_lup, _b);
}

void Task::task6_1()
{
	_ldlt = _matrix.LDLt();
}

void Task::task6_2()
{
	_ldlt_solution = Matrix::SLAE(_ldlt, _b);
}

void Task::task7(const double& w)
{
	_relax = _matrix.Relax(_b, w);
}
