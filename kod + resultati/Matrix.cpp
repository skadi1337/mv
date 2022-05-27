#include "Matrix.h"




size_t Matrix::main_element_by_column(const Matrix& matrix, size_t submatrix)
{
	size_t i = submatrix;

	double main_element = std::abs(matrix._data[i][i]);
	size_t main_index   = i;
	for (size_t k = i; k < matrix._rows; ++k)
	{
		double tmp = std::abs(matrix._data[k][i]);
		if (tmp > main_element)						// > COMPARISON
		{
			main_element = tmp;
			main_index   = k;
		}
	}

	return main_index;
}

//
//Complex::Matrix::Matrix(const ::Matrix& matrix)
//	: _rows		(matrix._rows)
//	, _columns	(matrix._columns)
//{ 
//	_data = std::vector<std::vector<std::complex<double>>>(_rows, std::vector<std::complex<double>>(_columns));
//
//	for (size_t i = 0; i < _rows; ++i)
//	{
//		for (size_t j = 0; j < _columns; ++j)
//		{
//			_data[i][j] = matrix._data[i][j];
//		}
//	}
//}

//Complex::Matrix::Matrix(size_t rows, size_t columns)
//	: _rows(rows)
//	, _columns(columns)
//{
//	_data = std::vector<std::vector<std::complex<double>>>(_rows, std::vector<std::complex<double>>(_columns));
//}
//
//Complex::Matrix Complex::Matrix::multiply(const Complex::Matrix & other) const
//{
//	auto result = Complex::Matrix(_rows, other._columns);
//
//	for (size_t i = 0; i < result._rows; i++)
//	{
//		for (size_t j = 0; j < result._columns; ++j)
//		{
//			std::complex<double> sum = 0;
//			for (size_t k = 0; k < _columns; ++k)
//			{
//				sum += _data[i][k] * other._data[k][j];
//			}
//
//			result._data[i][j] = sum;
//		}
//	}
//
//	return result;
//}
//
//Complex::Matrix Complex::Matrix::multiply(const std::complex<double>&value) const
//{
//	auto result = Complex::Matrix(_rows, _columns);
//
//	for (size_t i = 0; i < _rows; i++)
//	{
//		for (size_t j = 0; j < _columns; j++)
//		{
//			result._data[i][j] = _data[i][j] * value;
//		}
//	}
//
//	return result;
//}

//Complex::Matrix Complex::Matrix::substract(const Complex::Matrix& other) const
//{
//	auto result = Complex::Matrix(_rows, _columns);
//
//	for (size_t i = 0; i < _rows; ++i)
//	{
//		for (size_t j = 0; j < _columns; ++j)
//		{
//			result._data[i][j] = _data[i][j] - other._data[i][j];
//		}
//	}
//
//	return result;
//}

//double Complex::Matrix::norm() const
//{
//	double result = -std::numeric_limits<double>::infinity();
//	for (size_t i = 0; i < _rows; ++i)
//	{
//		double sum = 0;
//		for (size_t j = 0; j < _columns; ++j)
//		{
//			sum += std::abs(_data[i][j]);
//		}
//
//		result = std::max(result, sum);
//	}
//
//	return result;
//}
//
//size_t Complex::Matrix::rows() const
//{
//	return _rows;
//}
//
//size_t Complex::Matrix::columns() const
//{
//	return _columns;
//}


Matrix::Matrix(size_t rows, size_t columns)
	: _rows		(rows)
	, _columns  (columns)
{
	_data = std::vector<std::vector<matrix_data_t>>(_rows, std::vector<matrix_data_t>(_columns));
}

Matrix::Matrix(std::vector<std::vector<matrix_data_t>>&& data) noexcept
	: _rows		(data.size())
	, _columns  (data[0].size())
	, _data		(data)
{ }


Matrix::Matrix(size_t rows, size_t columns, const std::vector<std::vector<matrix_data_t>>& data)
	: _rows(rows)
	, _columns(columns)
	, _data(data)
{ }

Matrix::Matrix(size_t rows, size_t columns, std::vector<std::vector<matrix_data_t>>&& data) noexcept
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

Matrix Matrix::Generate(size_t rows, size_t columns, std::function<matrix_data_t()> generator)
{
	auto data = std::vector<std::vector<matrix_data_t>>(rows, std::vector<matrix_data_t>(columns));

	for (size_t i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < columns; j++)
		{
			data[i][j] = generator();
		}
	}

	return Matrix(rows, columns, data);
}

Matrix Matrix::FromFile(const std::filesystem::path& filepath)
{
	std::ifstream input(filepath);
	if (!input.is_open())
	{
		throw std::runtime_error("Could not open a file");
	}

	size_t rows, columns;
	input >> rows >> columns;

	auto data = std::vector<std::vector<matrix_data_t>>(rows, std::vector<matrix_data_t>(columns));

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

Matrix Matrix::transpose() const
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

Matrix Matrix::add(const Matrix& other) const
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

Matrix Matrix::substract(const Matrix& other) const
{
	auto result = this->copy();

	auto s = _rows - other._rows;

	for (size_t i = s; i < _rows; ++i)
	{
		for (size_t j = s; j < _columns; ++j)
		{
			result._data[i][j] = _data[i][j] - other._data[i - s][j - s];
		}
	}

	return result;
}

Matrix Matrix::multiply(const Matrix& other) const
{
	auto result = Matrix(_rows, other._columns);

	for (size_t i = 0; i < result._rows; i++)
	{
		for (size_t j = 0; j < result._columns; ++j)
		{
			matrix_data_t sum = 0;
			for (size_t k = 0; k < _columns; ++k)
			{
				sum += _data[i][k] * other._data[k][j];
			}

			result._data[i][j] = sum;
		}
	}

	return result;
}

Matrix Matrix::multiply(const matrix_data_t& value) const
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

Matrix Matrix::divide(const matrix_data_t& value) const
{
	auto result = Matrix(_rows, _columns);

	for (size_t i = 0; i < _rows; i++)
	{
		for (size_t j = 0; j < _columns; j++)
		{
			result._data[i][j] = _data[i][j] / value;
		}
	}

	return result;
}

//Complex::Matrix Complex::Matrix::SLAE(const Matrix& b) const
//{
//	auto matrix = this->copy();
//	auto solution = b.copy();
//
//
//	for (size_t i = 0; i < _rows - 1; ++i)
//	{
//		size_t main_index = main_element_by_column(matrix, i);
//
//		std::swap(matrix._data[i], matrix._data[main_index]);
//		std::swap(solution._data[i], solution._data[main_index]);
//
//
//		for (size_t k = i + 1; k < _rows; ++k)
//		{
//			matrix_data_t coef = matrix._data[k][i] / matrix._data[i][i];
//
//			for (size_t j = i + 1; j < _columns; ++j)
//			{
//				matrix._data[k][j] -= matrix._data[i][j] * coef;
//			}
//
//			solution._data[k][0] -= solution._data[i][0] * coef;
//
//		}
//
//	}
//
//
//	for (int i = _rows - 1; i >= 0; --i)
//	{
//		matrix_data_t sum = 0;
//		for (int j = _columns - 1; j > i; --j)
//		{
//			sum += matrix._data[i][j] * solution._data[j][0];
//		}
//
//		solution._data[i][0] = (solution._data[i][0] - sum) / matrix._data[i][i];
//
//	}
//
//	return solution;
//}

Matrix Matrix::SLAE(const Matrix& b) const
{
	auto matrix   = this->copy();
	auto solution = b.copy();


	for (size_t i = 0; i < _rows - 1; ++i)
	{
		size_t main_index = main_element_by_column(matrix,	i);

		std::swap(matrix._data[i], matrix._data[main_index]);
		std::swap(solution._data[i], solution._data[main_index]);


		for (size_t k = i + 1; k < _rows; ++k)
		{
			matrix_data_t coef = matrix._data[k][i] / matrix._data[i][i];

			for (size_t j = i + 1; j < _columns; ++j)
			{
				matrix._data[k][j] -= matrix._data[i][j] * coef;
			}

			solution._data[k][0] -= solution._data[i][0] * coef;

		}

	}


	for (int i = _rows - 1; i >= 0; --i)
	{
		matrix_data_t sum = 0;
		for (int j = _columns - 1; j > i; --j)
		{
			sum += matrix._data[i][j] * solution._data[j][0];
		}

		solution._data[i][0] = (solution._data[i][0] - sum) / matrix._data[i][i];
		
	}

	return solution;
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

matrix_data_t& Matrix::data(size_t row, size_t column)
{
	return this->_data[row][column];
}

Matrix& Matrix::operator=(Matrix&& other) noexcept
{
	_rows    = std::move(other._rows);
	_columns = std::move(other._columns);
	_data    = std::move(other._data);


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
