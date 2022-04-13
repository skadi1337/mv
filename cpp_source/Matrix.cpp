#include "Matrix.h"

std::vector<size_t> Matrix::swap_vector(size_t size)
{
	std::vector<size_t> swap(size);

	for (size_t i = 0; i < swap.size(); ++i)
	{
		swap[i] = i;
	}

	return swap;
}

size_t Matrix::main_element_by_column(const Matrix& matrix, size_t submatrix)
{
	size_t i = submatrix;

	double main_element = std::fabs(matrix._data[i][i]);
	size_t main_index   = i;
	for (size_t k = i; k < matrix._rows; ++k)
	{
		double tmp = std::fabs(matrix._data[k][i]);
		if (tmp > main_element)						// > COMPARISON
		{
			main_element = tmp;
			main_index   = k;
		}
	}

	return main_index;
}

size_t Matrix::main_element_by_row(const Matrix& matrix, size_t submatrix, const std::vector<size_t>& swap)
{
	size_t i = submatrix;

	double main_element = std::fabs(matrix._data[i][swap[i]]);
	size_t main_index	= i;
	for (size_t k = i; k < matrix._columns; ++k)
	{
		double tmp = std::fabs(matrix._data[i][swap[k]]);
		if (std::fabs(tmp - main_element) >= std::numeric_limits<double>::epsilon()) // epsilon COMPARISON
		{
			main_element = tmp;
			main_index	 = k;
		}
	}

	return main_index;
}

Matrix::Matrix(size_t rows, size_t columns)
	: _rows(rows)
	, _columns(columns)
{
	_data = std::vector<std::vector<double>>(rows, std::vector<double>(columns));
}

Matrix::Matrix(size_t rows, size_t columns, const std::vector<std::vector<double>>& data)
	: _rows(rows)
	, _columns(columns)
	, _data(data)
{ }

Matrix::Matrix(size_t rows, size_t columns, std::vector<std::vector<double>>&& data) noexcept
	: _rows(rows)
	, _columns(columns)
	, _data(std::move(data))
{ }

Matrix::Matrix(const Matrix& other)
	: _rows(other._rows)
	, _columns(other._columns)
	, _data(other._data)
{ }

Matrix::Matrix(Matrix&& other) noexcept
	: _rows(other._rows)
	, _columns(other._columns)
	, _data(std::move(other._data))
{ }

Matrix Matrix::Generate(size_t rows, size_t columns, std::function<double()> generator)
{
	auto data = std::vector<std::vector<double>>(rows, std::vector<double>(columns));

	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			data[i][j] = generator();
		}
	}

	return Matrix(rows, columns, data);
}

Matrix Matrix::FromFile(const std::string& filepath)
{
	std::ifstream input(filepath);
	if (!input.is_open())
	{
		throw std::runtime_error("Could not open a file");
	}

	size_t rows, columns;
	input >> rows >> columns;

	auto data = std::vector<std::vector<double>>(rows, std::vector<double>(columns));

	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			input >> data[i][j];
		}
	}

	input.close();


	return Matrix(rows, columns, data);
}

Matrix Matrix::FromSwap(const std::vector<size_t>& swap)
{
	auto data = std::vector<std::vector<double>>(swap.size(), std::vector<double>(swap.size(), 0));

	for (size_t i = 0; i < swap.size(); ++i)
	{
		data[i][swap[i]] = 1;
	}

	return Matrix(swap.size(), swap.size(), data);
}

Matrix Matrix::Transpose() const
{
	auto result = Matrix(_columns, _rows);

	for (size_t i = 0; i < _rows; ++i)
	{
		for (size_t j = 0; j < _columns; ++j)
		{
			result._data[j][i] = _data[i][j];
		}
	}

	return result;
}

Matrix Matrix::Add(const Matrix& other) const
{
	auto result = Matrix(_rows, _columns);

	for (size_t i = 0; i < _rows; ++i)
	{
		for (size_t j = 0; j < _columns; ++j)
		{
			result._data[i][j] = _data[i][j] + other._data[i][j];
		}
	}

	return result;
}

Matrix Matrix::Substract(const Matrix& other) const
{
	auto result = Matrix(_rows, _columns);

	for (size_t i = 0; i < _rows; ++i)
	{
		for (size_t j = 0; j < _columns; ++j)
		{
			result._data[i][j] = _data[i][j] - other._data[i][j];
		}
	}

	return result;
}

Matrix Matrix::Multiply(const Matrix& other) const
{
	auto result = Matrix(_rows, other._columns);

	for (size_t i = 0; i < result._rows; i++)
	{
		for (size_t j = 0; j < result._columns; ++j)
		{
			double sum = 0;
			for (size_t k = 0; k < _columns; ++k)
			{
				sum += _data[i][k] * other._data[k][j];
			}

			result._data[i][j] = sum;
		}
	}

	return result;
}

Matrix Matrix::Multiply(const double& value) const
{
	auto result = Matrix(_rows, _columns);

	for (size_t i = 0; i < _rows; i++)
	{
		for (size_t j = 0; j < _columns; j++)
		{
			result._data[i][j] = _data[i][j] * value;
		}
	}

	return result;
}

Matrix Matrix::Inverse() const
{
	auto matrix  = this->copy();
	auto inverse = IdentityMatrix(_rows);

	for (size_t i = 0; i < _rows; ++i)
	{
		double diag = matrix._data[i][i];
		for (size_t j = 0; j < _columns; ++j)
		{
			matrix._data[i][j] /= diag;
			inverse._data[i][j] /= diag;
		}

		for (size_t k = i + 1; k < _rows; ++k)
		{
			double coef = matrix._data[k][i];

			for (size_t j = 0; j < _columns; ++j)
			{
				matrix._data[k][j] -= matrix._data[i][j] * coef;
				inverse._data[k][j] -= inverse._data[i][j] * coef;
			}
		}
	}

	for (int i = _rows - 1; i >= 0; --i)
	{
		for (int j = i - 1; j >= 0; --j)
		{
			double coef = matrix._data[j][i];
			for (size_t k = 0; k < _columns; ++k)
			{
				inverse._data[j][k] -= inverse._data[i][k] * coef;
			}
		}
	}

	return inverse;
}

Matrix Matrix::SLAE(const Matrix& b) const
{
	auto matrix   = this->copy();
	auto solution = b.copy();


	for (size_t i = 0; i < _rows - 1; ++i)
	{
		size_t main_index = main_element_by_column(matrix, i);

		std::swap(matrix._data[i], matrix._data[main_index]);
		std::swap(solution._data[i], solution._data[main_index]);

		for (size_t k = i + 1; k < _rows; ++k)
		{
			double coef = matrix._data[k][i] / matrix._data[i][i];
			for (size_t j = i + 1; j < _columns; ++j)
			{
				matrix._data[k][j] -= matrix._data[i][j] * coef;
			}

			solution._data[k][0] -= solution._data[i][0] * coef;
		}

	}


	for (int i = _rows - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = _columns - 1; j > i; --j)
		{
			sum += matrix._data[i][j] * solution._data[j][0];
		}

		solution._data[i][0] = (solution._data[i][0] - sum) / matrix._data[i][i];
	}

	return solution;
}

Matrix Matrix::SLAE(const LUP_Decomposition& lup, const Matrix& b)
{
	auto solution = b.copy();

	for (size_t i = 0; i < lup.LU._rows; ++i)
	{
		double sum = 0;
		for (size_t j = 0; j < i; ++j)
		{
			sum += solution._data[j][0] * lup.LU._data[i][lup.P[j]];
		}
		solution._data[i][0] = solution._data[i][0] - sum;
	}


	for (int i = lup.LU._rows - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = lup.LU._columns - 1; j > i; --j)
		{
			sum += solution._data[j][0] * lup.LU._data[i][lup.P[j]];
		}
		solution._data[i][0] = (solution._data[i][0] - sum) / lup.LU._data[i][lup.P[i]];
	}


	auto res = Matrix(lup.P.size(), 1);
	for (size_t i = 0; i < lup.P.size(); ++i)
	{
		res._data[lup.P[i]][0] = solution._data[i][0];
	}


	return res;
}

Matrix Matrix::SLAE(const LDLt_Decomposition& ldlt, const Matrix& b)
{
	auto solution = b.copy();

	for (size_t i = 0; i < ldlt.L._rows; ++i)
	{
		double sum = 0;
		for (size_t j = 0; j < i; ++j)
		{
			sum += solution._data[j][0] * ldlt.L._data[i][j];
		}

		solution._data[i][0] = (solution._data[i][0] - sum) / ldlt.L._data[i][i];
	}
	
	for (size_t i = 0; i < ldlt.D.size(); ++i)
	{
		if (!ldlt.D[i])
		{
			solution._data[i][0] = -solution._data[i][0];
		}
	}

	for (int i = ldlt.L._rows - 1; i >= 0; --i)
	{
		double sum = 0;
		for (int j = ldlt.L._columns - 1; j > i; --j)
		{
			sum += solution._data[j][0] * ldlt.L._data[j][i];
		}
		solution._data[i][0] = (solution._data[i][0] - sum) / ldlt.L._data[i][i];
	}

	return solution;
}

LUP_Decomposition Matrix::LUP() const
{
	auto LU   = this->copy();
	auto swap = swap_vector(_columns);

	for (size_t i = 0; i < _rows - 1; ++i)
	{
		size_t swap_index = main_element_by_row(LU, i, swap);
		std::swap(swap[swap[i]], swap[swap_index]);

		for (size_t k = i + 1; k < _rows; ++k)
		{

			double coef = LU._data[k][swap[i]] / LU._data[i][swap[i]];
			for (size_t j = i + 1; j < _columns; ++j)
			{
				
				LU._data[k][swap[j]] -= coef * LU._data[i][swap[j]];
			}

			LU._data[k][swap[i]] = coef;
		}

	}


	return { LU, swap };
}

Matrix Matrix::U_from_LU(const Matrix& lu)
{
	auto l = lu.copy();

	for (size_t i = 0; i < lu._rows; ++i)
	{
		for (size_t j = 0; j < i; ++j)
		{
			l._data[i][j] = 0;
		}
	}

	return l;
}

Matrix Matrix::L_from_LU(const Matrix& lu)
{
	auto u = lu.copy();

	for (size_t i = 0; i < lu._rows; ++i)
	{
		u._data[i][i] = 1;
		for (size_t j = i + 1; j < lu._columns; ++j)
		{
			u._data[i][j] = 0;
		}
	}

	return u;
}

Matrix Matrix::Lt_from_LDLt(const Matrix& l)
{
	auto L = l.copy();

	for (size_t i = 0; i < l._rows; ++i)
	{
		for (size_t j = 0; j < i; ++j)
		{
			L._data[i][j] = 0;
		}

		for (size_t j = i + 1; j < l._columns; ++j)
		{
			L._data[i][j] = l._data[j][i];
		}
	}

	return L;
}

Matrix Matrix::D_from_LDLt(const std::vector<bool>& d)
{
	auto D = IdentityMatrix(d.size());

	for (size_t i = 0; i < d.size(); ++i)
	{
		D._data[i][i] = d[i] ? 1 : -1;
	}

	return D;
}

Matrix Matrix::L_from_LDLt(const Matrix& l)
{
	auto L = l.copy();

	for (size_t i = 0; i < l._rows; ++i)
	{
		for (size_t j = i + 1; j < l._columns; ++j)
		{
			L._data[i][j] = 0;	
		}
	}

	return L;
}

Matrix Matrix::ApplyRowSwap(const std::vector<size_t>& swap) const
{
	auto res = Matrix(_rows, _columns);
	for (size_t i = 0; i < _rows; ++i)
	{
		res._data[i] = _data[swap[i]];
	}

	return res;
}

Matrix Matrix::ApplyColumnSwap(const std::vector<size_t>& swap) const
{
	auto res = Matrix(_rows, _columns);
	for (size_t i = 0; i < _rows; ++i)
	{
		for (size_t j = 0; j < _columns; ++j)
		{
			res._data[i][j] = _data[i][swap[j]];
		}
	}

	return res;
}

LDLt_Decomposition Matrix::LDLt() const
{
	auto L = this->copy();
	auto D = std::vector<bool>(_rows);

	for (size_t i = 0; i < _rows; ++i)
	{
		double diag = std::sqrt(std::abs(L._data[i][i]));
		for (size_t k = i + 1; k < _rows; ++k)
		{
			double coef = L._data[i][k] / L._data[i][i];
			for (size_t j = i + 1; j < _columns; ++j)
			{
				L._data[k][j] -= L._data[i][j] * coef;
			}

			L._data[k][i] = coef * diag;
		}

		D[i] = L._data[i][i] >= 0;
		L._data[i][i] = diag;
	}

	return { L, D };
}

Relax_Result Matrix::Relax(const Matrix& _b, const double& w) const
{
	auto b		  =	_b.copy();
	auto matrix   = this->copy();
	auto solution = _b.copy();


	for (size_t i = 0; i < _rows; ++i)
	{
		double diag = matrix._data[i][i];
		b._data[i][0] /= matrix._data[i][i];
		for (size_t j = 0; j < _columns; ++j)
		{
			matrix._data[i][j] = -(matrix._data[i][j] / diag);
		}

		matrix._data[i][i] = 0;
		solution._data[i][0] = 1;
	}

	int num_of_iterations = 0;


	std::vector<double> norms;
	norms.reserve(40);


	double norm = std::numeric_limits<double>::infinity();
	do
	{
		++num_of_iterations;
		auto prev = solution.copy();

		for (size_t i = 0; i < _rows; ++i)
		{
			double sum = 0;
			for (size_t j = 0; j < _columns; ++j)
			{
				sum += matrix._data[i][j] * solution._data[j][0];
			}

			solution._data[i][0] = (1 - w) * solution._data[i][0] + w * (sum + b._data[i][0]);
		}

		norm = prev.Substract(solution).norm();
		norms.push_back(norm);

		if (num_of_iterations > MAX_RELAX_ITERATIONS)
		{
			return {norms, solution, -1};
		}

	} while (norm >= RELAX_LIMIT);


	return { norms, solution, num_of_iterations };
}

double Matrix::norm() const
{
	double result = -std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < _rows; ++i)
	{
		double sum = 0;
		for (size_t j = 0; j < _columns; ++j)
		{
			sum += std::abs(_data[i][j]);
		}

		result = std::max(result, sum);
	}

	return result;
}

Matrix Matrix::copy() const
{
	return Matrix(_rows, _columns, _data);
}

size_t Matrix::rows() const
{
	return _rows;
}

size_t Matrix::columns() const
{
	return _columns;
}

double& Matrix::data(size_t row, size_t column)
{
	return this->_data[row][column];
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
	_rows = std::move(other._rows);
	_columns = std::move(other._columns);
	_data = std::move(other._data);


	return *this;
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	for (size_t i = 0; i < m._rows; ++i)
	{
		for (size_t j = 0; j < m._columns; j++)
		{
			os << m._data[i][j] << " ";
		}

		os << "\n";
	}

	return os;
}

IdentityMatrix::IdentityMatrix(size_t size)
	: Matrix(size, size)
{
	for (size_t i = 0; i < size; ++i)
	{
		this->_data[i][i] = 1;
	}
}