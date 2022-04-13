#pragma once

#include <cmath>
#include <functional>

#include "Matrix.h"
#include "Random.h"


class Task
{
private:

	//default options

	static const double						MIN_RANDOM;
	static const double						MAX_RANDOM;

	static const std::function<double()>	RANDOM_GENERATOR;

	static constexpr size_t ROWS	= 256,
						    COLUMNS = ROWS;
	//---------------------------------------------------------------

public:

	size_t _rows,
		   _columns;

	Matrix	_matrix, // A
			_x,
			_b,
			_inverse,
			_gauss_solution,
			_lup_solution,
			_ldlt_solution;

	double _condition_number = 0;

	LUP_Decomposition	_lup;
	LDLt_Decomposition	_ldlt;
	Relax_Result		_relax;

	std::function<double()> _generator;

public:

	Task(size_t rows = ROWS, size_t columns = COLUMNS, std::function<double()> generator = RANDOM_GENERATOR);


	void task1();		// fill A

	void task2();		// fill x, b

	void task3_1();		// find inverse
	void task3_2();		// find conditional number 

	void task4();		// solve Gaus

	void task5_1();		// find	 LUP
	void task5_2();		// solve LUP

	void task6_1();		// find  LDLt
	void task6_2();		// solve LDLt

	void task7(const double& w = 0.75);	// solve Relax


};