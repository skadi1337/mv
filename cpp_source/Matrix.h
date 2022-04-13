#pragma once

#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <valarray>
#include <iostream>
#include <exception>
#include <functional>
#include <string_view>


class Matrix;
class IdentityMatrix;

struct LUP_Decomposition;
struct LDLt_Decomposition;
struct Relax_Result;


class Matrix
{

protected:

	static constexpr double RELAX_LIMIT			 = 1e-8;
	static constexpr int	MAX_RELAX_ITERATIONS = 3000;



	size_t				_rows;
	size_t				_columns;

	std::vector<std::vector<double>> _data;


	static std::vector<size_t> swap_vector(size_t size);
	
	static size_t main_element_by_column (const Matrix& matrix, size_t submatrix);
	static size_t main_element_by_row	 (const Matrix& matrix, size_t submatrix, const std::vector<size_t>& swap);
	

public:

	// ------------ Constructors --------------------------------------------------------------------------------------

	Matrix	(size_t rows = 0, size_t columns = 0);

	Matrix	(size_t rows, size_t columns, const std::vector<std::vector<double>>& data);
	Matrix	(size_t rows, size_t columns, std::vector<std::vector<double>>&& data) noexcept;

	Matrix	(const Matrix& other);
	Matrix	(Matrix&& other) noexcept;


	// ------------ Factories --------------------------------------------------------------------------------------

	static Matrix Generate	(size_t rows, size_t columns, std::function<double()> generator);
	static Matrix FromFile	(const std::string& filepath);
	static Matrix FromSwap	(const std::vector<size_t>& swap);

	// P shoud be already applied (e.g. ApplyColumnSwap to LU)
	static Matrix		U_from_LU	(const Matrix& lu),
						L_from_LU	(const Matrix& lu);

	static Matrix		L_from_LDLt  (const Matrix& l),
						Lt_from_LDLt (const Matrix& l),
						D_from_LDLt  (const std::vector<bool>& d);

	// ------------ Operations --------------------------------------------------------------------------------------

	Matrix Add			(const Matrix& other) const;
	Matrix Substract	(const Matrix& other) const;
	Matrix Multiply		(const Matrix& other) const;
	Matrix Multiply		(const double& value) const;


	Matrix Transpose () const;
	Matrix Inverse	 () const;
	

			Matrix SLAE	(const Matrix& b) const;
	static  Matrix SLAE	(const LUP_Decomposition&  lup, const Matrix& b);
	static  Matrix SLAE (const LDLt_Decomposition& ldlt, const Matrix& b);


	LUP_Decomposition   LUP			() const;


	Matrix ApplyRowSwap		(const std::vector<size_t>& swap) const;
	Matrix ApplyColumnSwap	(const std::vector<size_t>& swap) const;


	LDLt_Decomposition LDLt() const;



	Relax_Result Relax	(const Matrix& b, const double& w) const;


	double norm() const;
	Matrix copy() const;


	size_t  rows	() const;
	size_t	columns () const;
	double& data	(size_t row, size_t column);
	

	Matrix& operator = (Matrix&& other) noexcept;

	// ------------ Print --------------------------------------------------------------------------------------

	friend std::ostream&  operator << (std::ostream& os, const Matrix& m);

};


class IdentityMatrix : public Matrix
{
public:
	IdentityMatrix(size_t size);
};


struct LUP_Decomposition
{
	Matrix				LU;
	std::vector<size_t> P;
};

struct LDLt_Decomposition
{
	Matrix			   L;
	std::vector<bool>  D;
};

struct Relax_Result
{
	std::vector<double> norms;
	Matrix				solution;
	int					num_of_iterations;
};