#pragma once

#include <string>
#include <vector>
#include <utility>
#include <complex>
#include <fstream>
#include <iostream>
#include <exception>
#include <functional>
#include <filesystem>
#include <string_view>


class Matrix;
class IdentityMatrix;


using matrix_data_t = double;


class Matrix
{
public:

	
	Matrix	(size_t rows = 0, size_t columns = 0);
	Matrix	(std::vector<std::vector<matrix_data_t>>&& data) noexcept;

	Matrix	(size_t rows, size_t columns, const std::vector<std::vector<matrix_data_t>>& data);
	Matrix	(size_t rows, size_t columns, std::vector<std::vector<matrix_data_t>>&& data) noexcept;

	Matrix	(const Matrix& other);
	Matrix	(Matrix&& other) noexcept;


	

	static Matrix Generate	(size_t rows, size_t columns, std::function<matrix_data_t()> generator);
	static Matrix FromFile	(const std::filesystem::path& filepath);


	Matrix add			(const Matrix& other) const;
	Matrix substract	(const Matrix& other) const;
	Matrix multiply		(const Matrix& other) const;

	Matrix divide		(const matrix_data_t& value) const;
	Matrix multiply		(const matrix_data_t& value) const;


	Matrix transpose () const;

	
	Matrix SLAE	(const Matrix& b) const;
	


	double norm() const;
	Matrix copy() const;


	size_t  rows	() const;
	size_t	columns () const;
	matrix_data_t& data	(size_t row, size_t column);



	std::vector<std::vector<matrix_data_t>>&& movedata()
	{
		return std::move(_data);
	}
	

	Matrix& operator = (Matrix&& other) noexcept;

	friend std::ostream&  operator << (std::ostream& os, const Matrix& m);

	std::vector<std::vector<matrix_data_t>> _data;

protected:
	size_t				_rows;
	size_t				_columns;


	static size_t main_element_by_column (const Matrix& matrix, size_t submatrix);
};


class IdentityMatrix : public Matrix
{
public:
	IdentityMatrix(size_t size);
};


