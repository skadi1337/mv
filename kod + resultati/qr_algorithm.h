#pragma once

#include "Matrix.h"

constexpr double QR_ERROR = 1e-15;

void rotate_rows(Matrix& m, size_t i, size_t j, double cos, double sin)
{
	std::vector<matrix_data_t> new_row_i(m.columns()),
							   new_row_j(m.columns());

	for (size_t k = 0; k < m.columns(); ++k)
	{
		new_row_i[k] = cos * m._data[i][k] + sin * m._data[j][k];
		new_row_j[k] = cos * m._data[j][k] - sin * m._data[i][k];
 	}

	m._data[i] = std::move(new_row_i);
	m._data[j] = std::move(new_row_j);
}

void rotate_cols(Matrix& m, size_t i, size_t j, double cos, double sin)
{
	std::vector<matrix_data_t> new_col_i(m.rows()),
							   new_col_j(m.rows());

	for (size_t k = 0; k < m.rows(); ++k)
	{
		new_col_i[k] = cos * m._data[k][i] + sin * m._data[k][j];
		new_col_j[k] = cos * m._data[k][j] - sin * m._data[k][i];
	}

	for (size_t k = 0; k < m.rows(); ++k)
	{
		m._data[k][i] = new_col_i[k];
		m._data[k][j] = new_col_j[k];
	}
}


double square_norm(Matrix& m)
{
	double sum = 0;
	for (size_t i = 0; i < m.rows(); ++i)
	{
		sum += m.data(i, 0) * m.data(i, 0);
	}

	return std::sqrt(sum);
}

Matrix Hessenberg(const Matrix& A)
{
	auto m = A.copy();
	std::vector<std::vector<matrix_data_t>> d(m.rows() - 1, std::vector<matrix_data_t>(1));

	for (size_t j = 0; j < m.columns() - 2; ++j)
	{		
		for (size_t i = j + 1; i < m.rows(); ++i)
		{
			d[i - j - 1][0] = m.data(i, j);
		}

		auto x = Matrix(m.rows() - j - 1, 1, std::move(d));
		auto n = std::copysign(square_norm(x), -x.data(0, 0));
		if (std::abs(n) <= std::numeric_limits<double>::epsilon())
		{
			continue;
		}

		std::vector<std::vector<matrix_data_t>> w(m.rows() - j - 1, std::vector<matrix_data_t>(1));
		w[0][0] = n;

		auto ww = Matrix(std::move(w));
		auto v  = ww.substract(x);

		auto p = v.multiply(v.transpose()).divide(v.transpose().multiply(v).data(0, 0)).multiply(2);
		auto h = IdentityMatrix(m.rows()).substract(p);
	

		m = h.multiply(m).multiply(h);

		d = x.movedata();
	}

	return m;
}

struct angle
{
	double cos;
	double sin;
};


void QR(const Matrix& x, std::ostream& output = std::cout)
{
	auto m = Hessenberg(x);
	std::vector<std::complex<double>> evalues(m.rows());


	double current_error = 1;
	std::complex<double> root1, root2, discriminant;

	while (current_error > QR_ERROR)
	{
		current_error = -std::numeric_limits<double>::infinity();

		std::vector<angle> rotation_angles(m.rows());
		for (size_t i = 0; i < m.rows() - 1; ++i)
		{
			double divisor = std::sqrt(m.data(i, i) * m.data(i, i) + m.data(i + 1, i) * m.data(i + 1, i));
			if (std::abs(divisor) <= std::numeric_limits<double>::epsilon())
			{
				output << "Matrix has linearly dependent columns\n";
				break;
			}

			rotation_angles[i].cos = m.data(i, i)	  / divisor;
			rotation_angles[i].sin = m.data(i + 1, i) / divisor;

			rotate_rows(m, i, i + 1, rotation_angles[i].cos, rotation_angles[i].sin);
		}
		for (size_t i = 0; i < m.rows() - 1; ++i)
		{
			rotate_cols(m, i, i + 1, rotation_angles[i].cos, rotation_angles[i].sin);
		}


		bool is_complex = false;
		for (size_t i = 0; i < evalues.size(); i++)
		{
			if (is_complex)
			{
				double b = -(m.data(i - 1, i - 1) + m.data(i, i)),
					   c = -(m.data(i, i - 1) * m.data(i - 1, i)) + (m.data(i - 1, i - 1) * m.data(i, i));

				discriminant = b * b - 4 * c;

				root1 = (-b - std::sqrt(discriminant)) / (std::complex<double>)2;
				root2 = (-b + std::sqrt(discriminant)) / (std::complex<double>)2;

				current_error = std::max(std::abs(evalues[i] - root1), current_error);
				evalues[i] = root1;

				current_error = std::max(std::abs(evalues[i - 1] - root2), current_error);
				evalues[i - 1] = root2;

				is_complex = false;
				continue;
			}

			

			if (i == m.rows() - 1 || std::abs(m.data(i + 1, i)) <= std::numeric_limits<double>::epsilon())
			{
				root1		  = m.data(i, i);
				current_error = std::max(std::abs(evalues[i] - root1), current_error);

				evalues[i] = root1;
			}
			else
			{
				is_complex = true;
			}
			
		}

	}

	output << "EVALUES:\n";
	for (const auto& x : evalues)
	{
		output << x << "\n";
	}
}

