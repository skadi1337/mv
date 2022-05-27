#pragma once

#include "ComplexMatrix.h"


namespace power_method
{
	using complex::Matrix;
	using complex::matrix_data_t;


	constexpr double POWER_METHOD_ERROR			 = 1e-8;
	constexpr double POWER_METHOD_MAX_ITERATIONS = 400;

	Matrix evector1,
		   evector2;

	// [x1, x2] = [-x3]
	Matrix make_slae(Matrix& a, Matrix& b, Matrix& c) 
	{
		auto r = Matrix(a.rows(), 2);


		for (size_t i = 0; i < a.rows(); ++i)
		{
			r.data(i, 0) = a.data(i, 0);
			r.data(i, 1) = b.data(i, 0);

			c.data(i, 0) = -c.data(i, 0);
		}

		return r;
	}

	double findnorm(const Matrix& a, const Matrix& x, const std::complex<double>& evalue)
	{
		auto ax = a.multiply(x);

		auto ex = x.multiply(evalue);
		auto s  = ax.substract(ex);

		evector1 = x.copy();

		return s.norm();
	}

	double findnorm(const Matrix& a, const Matrix& x, const std::complex<double>& evalue1, const std::complex<double>& evalue2)
	{
		auto& ca = a;
		auto& cx = x;
		auto ax = ca.multiply(cx);

		{
			auto ex = cx.multiply(evalue1);
			auto s  = ax.substract(ex);

			evector2 = s.copy();
		}

		{
			auto ex = cx.multiply(evalue2);
			auto s  = ax.substract(ex);

			evector1 = s.copy();
		}


		auto nev1 = (ca.multiply(evector1)).substract((evector1.multiply(evalue1))),
			 nev2 = (ca.multiply(evector2)).substract((evector2.multiply(evalue2)));


		return std::max(nev1.norm(), nev2.norm());
	}

	Matrix init_vector(size_t size)
	{
		std::vector<std::vector<matrix_data_t>> data(size, std::vector<matrix_data_t>(1));
		for (auto& r : data)
			for (auto& c : r)
				c = 1;

		return Matrix(size, 1, std::move(data));
	}


	void base(const Matrix& A, std::ostream& output = std::cout, bool verbose = false)
	{
		auto x = init_vector(A.rows());

		auto x1 = x;
		auto x2 = A.multiply(x1);
		auto x3 = A.multiply(x2);

		std::complex<double> evalue1 = 0, 
							 evalue2 = 0;
		bool has_two_roots = false;


		double norm = findnorm(A, x3, evalue1);
		for (int i = 0; norm > POWER_METHOD_ERROR && i < POWER_METHOD_MAX_ITERATIONS; ++i)
		{
			if (verbose)
				output << "ITERATION #" << i << "\n";

			x1 = x3.divide(x3.norm());
			x2 = A.multiply(x1);
			x3 = A.multiply(x2);

			auto K   = make_slae(x1, x2, x3);
			auto K_T = K.transpose();

				 K  = K_T.multiply(K);
			auto b  = K_T.multiply(x3);

			auto c = K.SLAE(b);

			matrix_data_t discriminant = (c.data(1, 0) * c.data(1, 0) - (matrix_data_t)4 * c.data(0, 0));

			if (std::abs(discriminant) <= std::numeric_limits<double>::epsilon())
			{
				if (verbose)
					output << "ONE REAL ROOT (a1 = a2)  // this evalue might have two eigenvectors.";

				auto root1 = (-c.data(1, 0)) / (matrix_data_t)2;
				evalue1 = root1;

				norm = findnorm(A, x3, evalue1);

				has_two_roots = false;
			}
			else if (discriminant.real() > 0) 
			{
				auto root1 = (-c.data(1, 0) + std::sqrt(discriminant)) / std::complex<double>(2);
				auto root2 = (-c.data(1, 0) - std::sqrt(discriminant)) / std::complex<double>(2);

				if ((std::abs(root1) - std::abs(root2)) <= std::numeric_limits<double>::epsilon())
				{
					if (verbose)
						output << "TWO REAL ROOTS (|a1| = |a2|)";

					evalue1 = root1;
					evalue2 = root2;

					norm = findnorm(A, x3, evalue1, evalue2);
					has_two_roots = true;
				}
				else
				{
					if (verbose)
						output << "ONE REAL ROOT (|a1| > |a2|)";

					if (std::abs(root1) > std::abs(root2))
						evalue1 = root1;
					else
						evalue1 = root2;

					norm = findnorm(A, x3, evalue1);
					has_two_roots = false;
				}
			}
			else
			{
				if (verbose)
					output << "TWO COMPLEX ROOTS";

				auto dd = std::complex<double>(discriminant);

				std::complex<double> root1 = (-c.data(1, 0) + std::sqrt(dd)) / (std::complex<double>)2;
				std::complex<double> root2 = (-c.data(1, 0) - std::sqrt(dd)) / (std::complex<double>)2;

				evalue1 = root1;
				evalue2 = root2;

				norm = findnorm(A, x3, evalue1, evalue2);
				has_two_roots = true;
			}

			if (verbose)
			{
				output << "\n\nNORM = " << norm;
				output << "\n\nEVALUE  = " << evalue1;
				output << "\n\nEVECTOR = \n";
				output << evector1;


				if (has_two_roots)
				{
					output << "\n\nEVALUE  = " << evalue2;
					output << "\n\nEVECTOR = \n";
					output << evector2;
				}

				output << "\n\n";
			}
			
		}


		output << "\n\nNORM = " << norm;
		output << "\n\nEVALUE  = " << evalue1;
		output << "\n\nEVECTOR = \n";
		output << evector1;


		if (has_two_roots)
		{
			output << "\n\nEVALUE  = " << evalue2;
			output << "\n\nEVECTOR = \n";
			output << evector2;
		}

		output << "\n\n";

	}

}