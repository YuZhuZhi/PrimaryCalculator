#pragma once
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "DATA.h"
#define Max(a, b) (a > b) ? a : b;
#define Min(a, b) (a > b) ? b : a;

class Linear
{
public:
	class Vector
	{
	public:
		static void Input(std::vector<double>& x)
		{
			Data::Process::Input(x);
		}
		static void Input(std::vector<double>& x, int n)
		{
			Data::Process::Input(x, n);
		}
		static void Print(std::vector<double> x)
		{
			Data::Process::Print(x);
		}
		static double ScalarPro(std::vector<double> x, std::vector<double> y)
		{
			return (Data::Process::CroseSum(x, y));
		}
		static bool VectorPro(std::vector<double> x, std::vector<double> y, std::vector<double>& result)
		{
			if (x.size() == 4 && y.size() == 4) {
				result.clear();
				result.push_back(0);
				std::vector<std::vector<double> > temp;
				temp.resize(4, std::vector<double>(4));
				temp[2] = x, temp[3] = y;
				for (int i = 1; i <= 3; i++) result.push_back(Matrix::Cofactor(temp, 1, i));
				result[2] = -result[2];
				return true;
			}
			else return false;
		}
		static double Length(std::vector<double> x)
		{
			return (sqrt(Data::Process::Summary(x, 2)));
		}
		static double Angle(std::vector<double> x, std::vector<double> y)
		{
			return (acos((ScalarPro(x, y)) / (Length(x) * Length(y))));
		}

	};

	class Polynomial
	{
	public:
		static void Input(std::vector<double>& A, std::vector<double>& B)
		{
			Data::Process::Input(A);
			Data::Process::Input(B);
		}
		static void Multi(std::vector<double> A, std::vector<double> B, std::vector<double>& C)
		{
			int max = Max(A.size(), B.size());
			C.resize(max);

		}
	};

	class Equation
	{
	public:
		static void Print(std::vector<double>& result)
		{
			Data::Process::Print(result);
		}
		static void Print(std::vector<double>& result, int n)
		{
			for (int i = 1; i <= n; i++) {
				std::cout << "x_" << i << " = " << result[i] << std::endl;
			}
		}
		static bool Quadratic(double a, double b, double c, std::vector<double>& result)
		{
			result.clear();
			result.push_back(0);
			double Delta = b * b - 4 * a * c;
			if (Delta >= 0) {
				result.push_back(((-1) * b - sqrt(Delta)) / (2 * a));
				result.push_back(((-1) * b + sqrt(Delta)) / (2 * a));
				return true;
			}
			else return false;
		}
		static bool Cubic(double a, double b, double c, double d, std::vector<double>& result)
		{
			return true;
		}
		static bool Set(std::vector<std::vector<double> > coeff, std::vector<double> right, std::vector<double>& result)
		{
			result.clear();
			result.push_back(Matrix::Det(coeff));
			if (result[0] == 0) return false;
			else {
				for (int i = 1; i <= coeff.size() - 1; i++) {
					std::vector<std::vector<double> > temp = coeff;
					Matrix::Transpose(temp);
					temp[i] = right;
					result.push_back((Matrix::Det(temp)) / result[0]);
				}
				return true;
			}
		}
	};

	class Matrix
	{
	public:
		static void Input(std::vector<std::vector<double> >& matrix, int R, int C)
		{
			matrix.clear();
			double temp = 0;
			std::vector<double> vectp;
			for (int i = 1; i <= C + 1; i++) vectp.push_back(0);
			matrix.push_back(vectp);
			for (int r = 1; r <= R; r++) {
				vectp.resize(1, 0);
				for (int c = 1; c <= C; c++) {
					std::cin >> temp;
					vectp.push_back(temp);
				}
				matrix.push_back(vectp);
			}
		}
		static void Print(std::vector<std::vector<double> > matrix)
		{
			for (int r = 1; r <= matrix.size() - 1; r++) {
				for (int c = 1; c <= matrix[0].size() - 1; c++) {
					std::cout.width(6);
					std::cout.flags(std::ios::right);
					std::cout.precision(2);
					std::cout << std::fixed << matrix[r][c] << " ";
				}
				std::cout << std::endl << std::endl;
			}
		}
		static int Rank(std::vector<std::vector<double> > matrix)
		{
			return (matrix.size() - 1);
		}
		static int Row(std::vector<std::vector<double> > matrix)
		{
			return (matrix.size() - 1);
		}
		static int Column(std::vector<std::vector<double> > matrix)
		{
			return (matrix[1].size() - 1);
		}
		static void UnitMatrix(std::vector<std::vector<double> >& matrix, int n)
		{
			matrix.clear();
			matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) matrix[i][i] = 1;
		}
		static void UnitMatrix(std::vector<std::vector<double> >& matrix, int n, int num)
		{
			matrix.clear();
			matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) matrix[i][i] = num;
		}
		static void GausElmn(std::vector<std::vector<double> >& matrix)
		{
			int R = Row(matrix), C = Column(matrix);
			for (int c = 1; c <= C; c++) {
				int r = 1;
				for (r = c; r <= R; r++) {
					if (fabs(matrix[r][c]) >= 1e-8) break;
				}
				if (r <= R) {
					int temp = r;
					RowSwap(matrix, r, c);
					for (r = temp + 1; r <= R; r++) {
						if (matrix[r][c] == 0) {}
						else RowTrans(matrix, r, temp, -(matrix[r][c] / matrix[temp][c]));
					}
				}
			}
		}
		static void GausElmn(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2)
		{
			matrix2 = matrix1;
			int R = Row(matrix2), C = Column(matrix2);
			for (int c = 1; c <= C; c++) {
				int r = 1;
				for (r = c; r <= R; r++) {
					if (fabs(matrix2[r][c]) >= 1e-8) break;
				}
				if (r <= R) {
					int temp = r;
					RowSwap(matrix2, r, c);
					for (r = temp + 1; r <= R; r++) {
						if (matrix2[r][c] == 0) {}
						else RowTrans(matrix2, r, temp, -(matrix2[r][c] / matrix2[temp][c]));
					}
				}
			}
		}
		static double Det(std::vector<std::vector<double> > det)
		{
			int n = Rank(det);
			std::vector<std::vector<double> > temp = det;
			GausElmn(temp);
			double multi = 1;
			for (int i = 1; i <= n; i++) multi = multi * temp[i][i];
			return (multi);
		}
		static void Transpose(std::vector<std::vector<double> >& matrix)
		{
			int R = Column(matrix), C = Row(matrix);
			std::vector<std::vector<double> > temp;
			temp.resize(matrix[1].size(), std::vector<double>(matrix.size()));
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp[r][c] = matrix[c][r];
			}
			matrix = temp;
		}
		static void Transpose(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2)
		{
			int R = Column(matrix1), C = Row(matrix1);
			matrix2.clear();
			matrix2.resize(matrix1[1].size(), std::vector<double>(matrix1.size(), 0));
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) matrix2[r][c] = matrix1[c][r];
			}
		}
		static bool AddRight(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2)
		{
			if (matrix1.size() != matrix2.size()) return false;
			else {
				for (int r = 1; r <= Row(matrix1); r++) {
					std::vector<double> temp = matrix1[r];
					temp.insert(temp.end(), matrix2[r].begin() + 1, matrix2[r].end());
					matrix1[r] = temp;
				}
				return true;
			}
		}
		static bool AddRight(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3)
		{
			if (matrix1.size() != matrix2.size()) return false;
			else {
				matrix3.resize(matrix1.size(), std::vector<double>(matrix1[1].size() + matrix2[1].size() - 1));
				for (int r = 1; r <= Row(matrix1); r++) {
					std::vector<double> temp = matrix1[r];
					temp.insert(temp.end(), matrix2[r].begin() + 1, matrix2[r].end());
					matrix3[r] = temp;
				}
				return true;
			}
		}
		static bool AddDown(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2)
		{
			if (matrix1[1].size() != matrix2[1].size()) return false;
			else {
				for (int r = 1; r <= Row(matrix2); r++) matrix1.push_back(matrix2[r]);
				return true;
			}
		}
		static bool AddDown(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3)
		{
			if (matrix1[1].size() != matrix2[1].size()) return false;
			else {
				matrix3 = matrix1;
				for (int r = 1; r <= Row(matrix2); r++) matrix3.push_back(matrix2[r]);
				return true;
			}
		}
		static bool Plus(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2)
		{
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix1[r][c] = matrix1[r][c] + matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static bool Plus(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3)
		{
			matrix3.clear();
			matrix3.resize(matrix1.size(), std::vector<double>(matrix1[1].size(), 0));
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix3[r][c] = matrix1[r][c] + matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static bool Minus(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2)
		{
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix1[r][c] = matrix1[r][c] - matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static bool Minus(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3)
		{
			matrix3.clear();
			matrix3.resize(matrix1.size(), std::vector<double>(matrix1[1].size(), 0));
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix3[r][c] = matrix1[r][c] - matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static double Cofactor(std::vector<std::vector<double> > det, int R, int C)
		{
			std::vector<std::vector<double> > temp;
			temp.resize(det.size() - 1, std::vector<double>(det[1].size() - 1, 0));
			for (int r = 1; r <= Row(det); r++) {
				if (r < R) {
					for (int c = 1; c <= Column(det); c++) {
						if (c < C) temp[r][c] = det[r][c];
						if (c > C) temp[r][c - 1] = det[r][c];
					}
				}
				if (r > R) {
					for (int c = 1; c <= Column(det); c++) {
						if (c < C) temp[r - 1][c] = det[r][c];
						if (c > C) temp[r - 1][c - 1] = det[r][c];
					}
				}
			}
			return (Det(temp));
		}
		static bool Adjoint(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2)
		{
			matrix2 = matrix1;
			if (matrix1.size() != matrix1[1].size()) return false;
			else {
				for (int r = 1; r <= Rank(matrix1); r++) {
					for (int c = 1; c <= Rank(matrix1); c++) matrix2[r][c] = Cofactor(matrix1, c, r) * pow(-1, r + c);
				}
				return true;
			}
		}
		static bool Inverse(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2)
		{
			matrix2 = matrix1;
			double M = Det(matrix1);
			if (matrix1.size() != matrix1[1].size() || M == 0) return false;
			else {
				for (int r = 1; r <= Rank(matrix1); r++) {
					for (int c = 1; c <= Rank(matrix1); c++) matrix2[r][c] = (Cofactor(matrix1, c, r) / M) * pow(-1, r + c);
				}
				return true;
			}
		}
		static bool Multi(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3)
		{
			int R = Row(matrix1), C = Column(matrix2);
			matrix3.clear();
			matrix3.resize(matrix1.size(), std::vector<double>(matrix2[1].size(), 0));
			if (matrix1[1].size() != matrix2.size()) return false;
			else {
				Transpose(matrix2);
				for (int r = 1; r <= R; r++) {
					for (int c = 1; c <= C; c++) matrix3[r][c] = Vector::ScalarPro(matrix1[r], matrix2[c]);
				}
				return true;
			}
		}
		static bool Power(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2, int power)
		{
			if (matrix1.size() != matrix1[1].size()) return false;
			else {
				matrix2 = matrix1;
				if (power > 1) {
					for (int i = 2; i <= power; i++) {
						Multi(matrix2, matrix1, matrix2);
					}
				}
				if (power == 1) return true;
				if (power == 0) {
					UnitMatrix(matrix2, matrix1.size() - 1);
				}
				if (power < 0) {
					Inverse(matrix1, matrix2);
					std::vector<std::vector<double> > matrix3 = matrix2;
					for (int i = 2; i <= -power; i++) {
						Multi(matrix2, matrix3, matrix2);
					}
				}
				return true;
			}
		}
		static bool Eigen()
		{

		}

	protected:
		static void RowTrans(std::vector<std::vector<double> >& det, int r, double k)
		{
			int n = det[1].size() - 1;
			for (int i = 1; i <= n; i++) det[r][i] = k * det[r][i];
		}
		static void RowTrans(std::vector<std::vector<double> >& det, int r1, int r2, double k)
		{
			int n = det[1].size() - 1;
			for (int i = 1; i <= n; i++) det[r1][i] = det[r1][i] + k * det[r2][i];
		}
		static void RowSwap(std::vector<std::vector<double> >& det, int r1, int r2)
		{
			det[r1].swap(det[r2]);
		}
	};

};