#include "Task.h"
#include "Timer.h"

#include <algorithm>

void example ()
{
	Task task(3, 3, default_int_generator(-5, 5));

	std::cout << "Base ==================================\n";
	task._matrix = Matrix::FromFile("tests\\example.txt");
	task.task2();
	task.task3_1();
	task.task3_2();
	task.task4();
	task.task5_1();
	task.task5_2();
	std::cout << "----------------------\nA:\n" << task._matrix << "----------------------\nx:\n" << task._x << "----------------------\nb:\n" << task._b;
	std::cout << "----------------------\nInverse:\n" << task._inverse;
	std::cout << "----------------------\nConditional number: " << task._condition_number << "\n";
	std::cout << "----------------------\nSLAE:\n" << task._gauss_solution << "\n";

	task._lup.LU = task._lup.LU.ApplyColumnSwap(task._lup.P);
	auto p = Matrix::FromSwap(task._lup.P);
	auto l = Matrix::L_from_LU(task._lup.LU);
	auto u = Matrix::U_from_LU(task._lup.LU);

	auto lu = l.Multiply(u);
	auto lup = lu.Multiply(p);
	std::cout << "----------------------\nL:\n" << l;
	std::cout << "----------------------\nU:\n" << u;
	std::cout << "----------------------\nP:\n" << p;
	std::cout << "----------------------\nLUP\n" << lup;
	std::cout << "----------------------\nLUP Solution:\n" << task._lup_solution;


	std::cout << "RELAX ========================================\n";

	task._matrix = Matrix::FromFile("tests\\example_relax.txt");
	std::vector<std::vector<double>> b = { { 1}, {-3}, {0} };
	task._x = Matrix(3, 1, b);
	task._b = task._matrix.Multiply(task._x);

	std::cout << "----------------------\nA:\n" << task._matrix << "----------------------\nx:\n" << task._x << "----------------------\nb:\n" << task._b;

	task.task7();
	std::cout << "\n\n";
	std::cout << "Relax solution:\n" << task._relax.solution << "\n\n";
	std::cout << "iterations: " << task._relax.num_of_iterations << "\n";

	std::cout << "LDLT=====================================================================\n";
	task._matrix = Matrix::FromFile("tests\\example_ldlt.txt");
	task.task2();
	std::cout << "----------------------\nA:\n" << task._matrix << "----------------------\nx:\n" << task._x << "----------------------\nb:\n" << task._b;
	auto ldlt = task._matrix.LDLt();

	auto ll = Matrix::L_from_LDLt(ldlt.L);
	auto d = Matrix::D_from_LDLt(ldlt.D);
	auto llt = Matrix::Lt_from_LDLt(ldlt.L);

	std::cout << "----------------------\nL:\n" << ll;
	std::cout << "----------------------\nD:\n" << d;
	std::cout << "----------------------\nLt:\n" << llt; 
	std::cout << "----------------------\nLDLt:\n" << ll.Multiply(d).Multiply(llt);
	std::cout << "----------------------\nLDLt SLAE:\n" << Matrix::SLAE(ldlt, task._b);
}

void experiment (int num_of_experiments = 100)
{
	Timer timer;
	std::ofstream output("results\\otchet.txt");

	// 8.1
	double sum_of_conditional_numbers = 0;

	double avg_conditional_number = 0,
		   max_conditional_number = -std::numeric_limits<double>::infinity(),
		   min_conditional_number =  std::numeric_limits<double>::infinity();

	Matrix max_conditional_matrix;


	//8.2
	double avg_time_inverse   = 0,	
		   total_time_inverse = 0;

	//8.3 
	double	min_norm_gaus =  std::numeric_limits<double>::infinity(),
			max_norm_gaus = -std::numeric_limits<double>::infinity(),
			avg_norm_gaus = 0,
			sum_norm_gaus = 0,

			min_norm_lup =  std::numeric_limits<double>::infinity(),
			max_norm_lup = -std::numeric_limits<double>::infinity(),
			avg_norm_lup = 0,
			sum_norm_lup = 0,

			min_norm_ldlt =  std::numeric_limits<double>::infinity(),
			max_norm_ldlt = -std::numeric_limits<double>::infinity(),
			avg_norm_ldlt = 0,
			sum_norm_ldlt = 0,

			min_norm_relax =  std::numeric_limits<double>::infinity(),
			max_norm_relax = -std::numeric_limits<double>::infinity(),
			avg_norm_relax = 0,
			sum_norm_relax = 0;

	double norm;

	//8.4
	double avg_time_gaus	= 0,
		   total_time_gaus	= 0;

	//8.5
	double avg_time_lup_make   = 0,
		   total_time_lup_make = 0;

	//8.6
	double avg_time_lup_solve	= 0,
		   total_time_lup_solve = 0;

	//8.7
	double avg_time_ldlt_solve	 = 0,
		   total_time_ldlt_solve = 0;

	//8.8
	double avg_time_relax_solve	  = 0,
		   total_time_relax_solve = 0;

	//8.9
	int min_it_cnt = std::numeric_limits<int>::max(),
		max_it_cnt = std::numeric_limits<int>::min(),
		avg_it_cnt = 0,
		sum_it_cnt = 0;

	std::vector<double> errors;

	for (int i = 0; i < num_of_experiments; ++i)
	{
		Task task;
		task.task1();
		task.task2();

		timer.start();
		task.task3_1();
		timer.stop();

		total_time_inverse += timer.elapsed();

		task.task3_2();
		sum_of_conditional_numbers += task._condition_number;
		min_conditional_number = std::min(min_conditional_number, task._condition_number);
		if (task._condition_number > max_conditional_number)
		{
			max_conditional_number = task._condition_number;
			max_conditional_matrix = task._matrix.copy();
		}


		// Gaus
		timer.start();
		task.task4();
		timer.stop();

		total_time_gaus += timer.elapsed();

		norm = task._gauss_solution.Substract(task._x).norm();

		min_norm_gaus = std::min(min_norm_gaus, norm);
		max_norm_gaus = std::max(max_norm_gaus, norm);
		sum_norm_gaus += norm;

		// LUP
		timer.start();
		task.task5_1();
		timer.stop();

		total_time_lup_make += timer.elapsed();


		timer.start();
		task.task5_2();
		timer.stop();

		total_time_lup_solve += timer.elapsed();

		norm = task._lup_solution.Substract(task._x).norm();

		min_norm_lup = std::min(min_norm_lup, norm);
		max_norm_lup = std::max(max_norm_lup, norm);
		sum_norm_lup += norm;

		// LDLt

		task.task6_1();

		timer.start();
		task.task6_2();
		timer.stop();

		total_time_ldlt_solve += timer.elapsed();

		norm = task._ldlt_solution.Substract(task._x).norm();

		min_norm_ldlt = std::min(min_norm_ldlt, norm);
		max_norm_ldlt = std::max(max_norm_ldlt, norm);
		sum_norm_ldlt += norm;

		// Relax

		timer.start();
		task.task7();
		timer.stop();

		total_time_relax_solve += timer.elapsed();

		norm = task._relax.solution.Substract(task._x).norm();

		min_norm_relax = std::min(min_norm_relax, norm);
		max_norm_relax = std::max(max_norm_relax, norm);
		sum_norm_relax += norm;

		min_it_cnt = std::min(min_it_cnt, task._relax.num_of_iterations);
		max_it_cnt = std::max(max_it_cnt, task._relax.num_of_iterations);
		sum_it_cnt += task._relax.num_of_iterations;

		if (i == num_of_experiments - 1)
		{
			errors.assign(task._relax.norms.begin(), task._relax.norms.end());
		}
		
	}


	avg_time_gaus		   = total_time_gaus			/ num_of_experiments;
	avg_time_inverse	   = total_time_inverse			/ num_of_experiments;
	avg_conditional_number = sum_of_conditional_numbers / num_of_experiments;

	avg_time_lup_make    = total_time_lup_make    / num_of_experiments;
	avg_time_lup_solve   = total_time_lup_solve   / num_of_experiments;
	avg_time_ldlt_solve  = total_time_ldlt_solve  / num_of_experiments;
	avg_time_relax_solve = total_time_relax_solve / num_of_experiments;

	avg_norm_gaus  = sum_norm_gaus  / num_of_experiments;
	avg_norm_lup   = sum_norm_lup   / num_of_experiments;
	avg_norm_ldlt  = sum_norm_ldlt  / num_of_experiments;
	avg_norm_relax = sum_norm_relax / num_of_experiments;

	avg_it_cnt	   = sum_it_cnt		/ num_of_experiments;


	output << "\n min conditional_number: " << min_conditional_number;
	output << "\n max conditional_number: " << max_conditional_number;
	output << "\n avg conditional_number: " << avg_conditional_number;
	
	output << "\n\n avg time gaus:\t\t"    << avg_time_gaus << " ms";
	output << "\n avg time inverse:\t"     << avg_time_inverse << " ms";
	output << "\n avg time lup_make:\t"    << avg_time_lup_make << " ms";
	output << "\n avg time lup_solve:\t"   << avg_time_lup_solve << " ms";
	output << "\n avg time ldlt_solve:\t"  << avg_time_ldlt_solve << " ms";
	output << "\n avg time relax_solve:\t" << avg_time_relax_solve << " ms";
	
	output << "\n\n min norm gaus: " << min_norm_gaus;
	output << "\n max norm gaus: "   << max_norm_gaus;
	output << "\n avg norm gaus: "   << avg_norm_gaus;
	
	output << "\n\n min norm lup: " << min_norm_lup;
	output << "\n max norm lup: "   << max_norm_lup;
	output << "\n avg norm lup: "   << avg_norm_lup;
	
	output << "\n\n min norm ldlt: " << min_norm_ldlt;
	output << "\n max norm ldlt: "   << max_norm_ldlt;
	output << "\n avg norm ldlt: "   << avg_norm_ldlt;
	
	output << "\n\n min norm relax: " << min_norm_relax;
	output << "\n max norm relax: "   << max_norm_relax;
	output << "\n avg norm relax: "   << avg_norm_relax;
	
	output << "\n\n min relax iterations: "	<< min_it_cnt;
	output << "\n max relax iterations: "	<< max_it_cnt;
	output << "\n avg relax iterations: "	<< avg_it_cnt;

	output.close();

	output.open("results\\max_conditional_matrix.txt");

	output << max_conditional_matrix.rows() << " " << max_conditional_matrix.columns() << "\n";
	output << max_conditional_matrix;

	output.close();
}

constexpr double N = 10;

const std::vector<std::vector<double>> data1 = { {N * N + 15, N - 1		 , -1		, -2		},
												 {N - 1		, -15 - N * N, -N + 4	, -4		},
												 {-1		, -N + 4	 , N * N + 8, -N		},
												 {-2		, -4		 , -N		, N * N + 10} };

const std::vector<std::vector<double>> data2 = { {1, 1 + N, 2 + N, 3 + N, 4 + N, 5 + N, 6 + N, 7 + N},
												  {100 * N, 1000 * N, 10000 * N, 100000 * N, -1000 * N, -10000 * N, -100000 * N, 1},
												  {N, -1 + N, -2 + N, -3 + N, -4 + N, -5 + N, -6 + N, -7 + N},
												  {N - 1000, 10 * N - 1000, 100 * N - 1000, 1000 * N - 1000, 10000 * N - 1000, -N, -N + 1, -N + 2},
												  {-2 * N, 0, -1, -2,-3,-4,-5,-6},
												  {N - 2019, -N + 2020, N-2021,-N+2022,N-2023,-N+2024,N-2025,-N+2026},
												  {2*N-2000,4*N-2005,8*N-2010,16*N-2015,32*N-2020,2019*N,-2020*N,2021*N},
												  {1020-2*N,-2924+896*N,1212+9808*N,-2736+98918*N,1404-11068*N,-1523-8078*N,2625-102119*N,-1327+1924*N} };

void task9(const Matrix& m, const std::string& file)
{
	Task task(m.rows(), m.columns());
	task._matrix = m.copy();

	task.task2();
	task.task3_1();
	task.task3_2();
	task.task4();
	task.task5_1();
	task.task5_2();
	task.task6_1();
	task.task6_2();
	task.task7();

	std::ofstream output(file);


	output << "----------------------\nConditional number: " << task._condition_number << "\n";
	output << "\n----------------------\nSolutions norms:\n";

	output << "Gaus norm:  " <<  task._gauss_solution.Substract(task._x).norm();
	output << "\nLUP norm:   " << task._lup_solution.Substract(task._x).norm();
	output << "\nLDLt norm:  " << task._ldlt_solution.Substract(task._x).norm();
	output << "\nRelax norm: " << task._relax.solution.Substract(task._x).norm();


	output << "\n\n----------------------\nA:\n" << task._matrix << "----------------------\nx:\n" << task._x << "----------------------\nb:\n" << task._b;
	output << "----------------------\nInverse:\n" << task._inverse;
	output << "----------------------\nConditional number: " << task._condition_number << "\n";
	output << "----------------------\nGauss solution:\n" << task._gauss_solution << "\n";

	auto lupp = task._lup.LU.ApplyColumnSwap(task._lup.P);
	auto p = Matrix::FromSwap(task._lup.P);
	auto l = Matrix::L_from_LU(lupp);
	auto u = Matrix::U_from_LU(lupp);

	auto lu = l.Multiply(u);
	auto lup = lu.Multiply(p);
	output << "----------------------\nL:\n" << l;
	output << "----------------------\nU:\n" << u;
	output << "----------------------\nP:\n" << p;
	output << "----------------------\nLUP\n" << lup;
	output << "----------------------\nLUP Solution:\n" << task._lup_solution;
	output << "----------------------\nRelax Solution:\n" << task._relax.solution << "\n\n";
	output << "iterations: " << task._relax.num_of_iterations << "\n";
	auto ldlt = task._matrix.LDLt();

	auto ll = Matrix::L_from_LDLt(ldlt.L);
	auto d = Matrix::D_from_LDLt(ldlt.D);
	auto llt = Matrix::Lt_from_LDLt(ldlt.L);

	output << "----------------------\nL:\n" << ll;
	output << "----------------------\nD:\n" << d;
	output << "----------------------\nLt:\n" << llt;
	output << "----------------------\nLDLt:\n" << ll.Multiply(d).Multiply(llt);
	output << "----------------------\nLDLt Solution:\n" << Matrix::SLAE(ldlt, task._b);




	output.close();
}

void task10_1(const Matrix& m, const std::string& file, int num_of_experiments = 10, const double& delta = 1e-5)
{
	std::ofstream output(file);

	Task task(m.rows(), m.columns());
	task._matrix = m.copy();
	task.task2();

	for (int i = 0; i < num_of_experiments; ++i)
	{
		for (size_t j = 0; j < task._b.rows(); ++j)
		{
			task._b.data(j, 0) += delta;
		}

		auto gaus   = task._matrix.SLAE(task._b);
		double norm = task._x.Substract(gaus).norm();

		output << norm << "\t\t\tdelta=" << std::to_string((i+1ll) * delta) << "\n";
	}

	output.close();
}

void task10_2(const Matrix& m, const std::string& file)
{
	std::ofstream output(file);

	Task task(m.rows(), m.columns());
	task._matrix = m.copy();
	task.task2();

	std::vector<double> params = { 0.8, 1.0, 1.2 };

	output << "\n\n=================================================\n\n";

	for (size_t i = 0; i < params.size(); ++i)
	{
		task.task7(params[i]);

		output << "w = " << std::to_string(params[i]) << " {";
		for (size_t i = 0; i < task._relax.norms.size() && i < 500; ++i)
		{
			output << task._relax.norms[i] << ", ";
		}
		output << "}\n\n";
	}

	output.close();
}


int main()
{
	Matrix a1(4, 4, data1);
	Matrix a2(8, 8, data2);

	a2 = a2.Transpose().Multiply(a2);

	std::cout << "Preparing otchet ...\n";
	experiment(100);
	

	auto m = Matrix::FromFile("results\\max_conditional_matrix.txt");

	task9(a1, "results\\A1.txt");
	task9(a2, "results\\A2.txt");

	task10_1(a2, "results\\research_A2.txt");
	task10_1(m, "results\\research_max_conditional.txt");

	task10_2(a2, "results\\errors_A2.txt");
	task10_2(m, "results\\errors_max_conditional.txt");

	std::cout << "Otchet done\n";

	//example();

	return 0;
}