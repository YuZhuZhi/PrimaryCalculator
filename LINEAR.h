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
	class Vector //��������
	{
	public:
		static void Input(std::vector<double>& x) //��������
		{
			Data::Process::Input(x);
		}
		static void Input(std::vector<double>& x, int n) //ָ��ά������������
		{
			Data::Process::Input(x, n);
		}
		static void Print(std::vector<double> x) //�������
		{
			std::cout << "(";
			for (std::vector<double>::iterator vdit = x.begin() + 1; vdit != x.end() - 1; vdit++) {
				std::cout << *vdit << ", ";
			}
			std::cout << *(x.end() - 1) << ")" << std::endl << std::endl;
		}
		static double ScalarPro(std::vector<double> x, std::vector<double> y) //����������������
		{
			return (Data::Process::CroseSum(x, y)); //���������ݼ��Ľ����
		}
		static bool VectorPro(std::vector<double> x, std::vector<double> y, std::vector<double>& result) //������������������result������������
		{
			if (x.size() == 4 && y.size() == 4) {
				result.clear();
				result.push_back(0);
				std::vector<std::vector<double> > temp;
				temp.resize(4, std::vector<double>(4)); //����һ����������ʽ
				temp[2] = x, temp[3] = y;
				for (int i = 1; i <= 3; i++) result.push_back(Matrix::Cofactor(temp, 1, i)); //�����Ӧ����ʽ
				result[2] = -result[2]; //�ڶ�ά��ֵȡ�෴��
				return true;
			}
			else return false;
		}
		static double Length(std::vector<double> x) //���������ĳ���
		{
			return (sqrt(Data::Process::Summary(x, 2))); //���������ݼ���ƽ����
		}
		static double Angle(std::vector<double> x, std::vector<double> y) //������������н�
		{
			return (acos((ScalarPro(x, y)) / (Length(x) * Length(y))));
		}

	};

	class Polynomial //����ʽ
	{
	public:
		static void Input(std::vector<double>& A, std::vector<double>& B) //����ʽϵ������
		{
			Data::Process::Input(A);
			Data::Process::Input(B);
		}
		static void Multi(std::vector<double> A, std::vector<double> B, std::vector<double>& C) //������ʽ���
		{
			int max = Max(A.size(), B.size());
			C.resize(max);

		}
	};

	class Equation //�ⷽ��
	{
	public:
		static void Print(std::vector<double> result) //������̵Ķ����
		{
			int n = Data::Process::Capacity(result);
			for (int i = 1; i <= n; i++) {
				std::cout << "x_"  << i << " = " << result[i] << std::endl;
			}
		}
		static bool Quadratic(double a, double b, double c, std::vector<double>& result) //����η��̣���result��װ��
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
		static bool Cubic(double a, double b, double c, double d, std::vector<double>& result) //�����η��̣���result��װ��
		{
			return true;
		}
		static bool Set(std::vector<std::vector<double> > coeff, std::vector<double> right, std::vector<double>& result) //�ⷽ���飬����result��װ��
		{
			result.clear();
			result.push_back(Matrix::Det(coeff)); //0λ������ϵ������ʽ
			if (result[0] == 0) return false;
			else {
				for (int i = 1; i <= Matrix::Rank(coeff); i++) {
					std::vector<std::vector<double> > temp = coeff;
					Matrix::Transpose(temp);
					temp[i] = right;
					result.push_back((Matrix::Det(temp)) / result[0]);
				}
				return true;
			}
		}
	};

	class Matrix //����ʽ���������
	{
	public:
		static void Input(std::vector<std::vector<double> >& matrix, int R, int C) //��������R������C�ľ�������
		{
			matrix.clear();
			double temp = 0;
			std::vector<double> vectp;
			for (int i = 1; i <= C + 1; i++) vectp.push_back(0);
			matrix.push_back(vectp); //��άvector��0������0Ԫ��
			for (int r = 1; r <= R; r++) {
				vectp.resize(1, 0);
				for (int c = 1; c <= C; c++) {
					std::cin >> temp;
					vectp.push_back(temp);
				}
				matrix.push_back(vectp);
			}
		}
		static void Print(std::vector<std::vector<double> > matrix) //��������
		{
			for (int r = 1; r <= Row(matrix); r++) {
				for (int c = 1; c <= Column(matrix); c++) {
					std::cout.width(6); //�ܼ���Ԫ�ؽ�С��ÿ��Ԫ����6��Ŀռ����
					std::cout.flags(std::ios::right); //�̶������λС�������Ҷ���
					std::cout.precision(2); //�̶������λС��
					std::cout << std::fixed << matrix[r][c] << " ";
				}
				std::cout << std::endl << std::endl;
			}
		}
		static int Rank(std::vector<std::vector<double> > matrix) //���ط���Ľ���
		{
			return (matrix.size() - 1);
		}
		static int Row(std::vector<std::vector<double> > matrix) //���ؾ��������
		{
			return (matrix.size() - 1);
		}
		static int Column(std::vector<std::vector<double> > matrix) //���ؾ��������
		{
			return (matrix[1].size() - 1);
		}
		static void UnitMatrix(std::vector<std::vector<double> >& matrix, int n) //���ɵ�λ����
		{
			matrix.clear();
			matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) matrix[i][i] = 1;
		}
		static void UnitMatrix(std::vector<std::vector<double> >& matrix, int n, int num) //���ɵ�λ�����num��
		{
			matrix.clear();
			matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) matrix[i][i] = num;
		}
		static void GausElmn(std::vector<std::vector<double> >& matrix) //��ԭ�����˹��Ԫ�������н�����
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
		static void GausElmn(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2) //��ԭ�����˹��Ԫ�����ɵ��н����δ���matrix2
		{
			matrix2 = matrix1;
			GausElmn(matrix2);
		}
		static void RowSim(std::vector<std::vector<double> >& matrix) //��ԭ���������������
		{
			int R = Row(matrix), C = Column(matrix);
			GausElmn(matrix);
			for (int r = 1; r <= R; r++) {
				int c = 1;
				for (c = 1; c <= C; c++) {
					if (fabs(matrix[r][c]) >= 1e-6) break;
				}
				RowTrans(matrix, r, 1.0 / matrix[r][c]);
			}
			for (int r = R; r >= 1; r--) {
				int c = C;
				for (c = 1; c <= C; c++) {
					if (fabs(matrix[r][c]) >= 1e-6) break;
				}
				for (int r2 = r - 1; r2 >= 1; r2--) {
					RowTrans(matrix, r2, r, -matrix[r2][c]);
				}
			}
		}
		static void RowSim(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2) //��ԭ�������ɵ�������δ���matrix2
		{
			matrix2 = matrix1;
			RowSim(matrix2);
		}
		static double Det(std::vector<std::vector<double> > det) //��������ʽ��ֵ
		{
			int n = Rank(det);
			std::vector<std::vector<double> > temp = det;
			GausElmn(temp); //��˹��ԪΪ�Ͻ�����
			double multi = 1;
			for (int i = 1; i <= n; i++) multi = multi * temp[i][i];
			return (multi);
		}
		static void Transpose(std::vector<std::vector<double> >& matrix) //��ԭ����ת��
		{
			int R = Column(matrix), C = Row(matrix);
			std::vector<std::vector<double> > temp;
			temp.resize(matrix[1].size(), std::vector<double>(matrix.size()));
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp[r][c] = matrix[c][r];
			}
			matrix = temp;
		}
		static void Transpose(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2) //��ԭ����ת�ú����matrix2
		{
			matrix2 = matrix1;
			Transpose(matrix2);
		}
		static bool AddRight(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2) //��ԭ�����ұ߸���ͬ��������
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
		static bool AddRight(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3) //��ԭ�����ұ߸���ͬ������������matrix3
		{
			matrix3 = matrix1;
			AddRight(matrix3, matrix2);
		}
		static bool AddDown(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2) //��ԭ�������渽��ͬ��������
		{
			if (matrix1[1].size() != matrix2[1].size()) return false;
			else {
				for (int r = 1; r <= Row(matrix2); r++) matrix1.push_back(matrix2[r]);
				return true;
			}
		}
		static bool AddDown(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3) //��ԭ�������渽��ͬ������������matrix3
		{
			matrix3 = matrix1;
			AddDown(matrix3, matrix2);
		}
		static bool Plus(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2) //��matrix1����matirx2
		{
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix1[r][c] = matrix1[r][c] + matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static bool Plus(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3) //��matrix1����matirx2�Ľ������matirx3
		{
			matrix3 = matrix1;
			Plus(matrix3, matrix2);
		}
		static bool Minus(std::vector<std::vector<double> >& matrix1, std::vector<std::vector<double> > matrix2) //��matrix1��ȥmatirx2
		{
			if (matrix1.size() == matrix2.size() && matrix1[1].size() == matrix2[1].size()) {
				for (int r = 1; r <= Row(matrix1); r++) {
					for (int c = 1; c <= Column(matrix1); c++) matrix1[r][c] = matrix1[r][c] - matrix2[r][c];
				}
				return true;
			}
			else return false;
		}
		static bool Minus(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3) //��matrix1��ȥmatirx2�Ľ������matirx3
		{
			matrix3 = matrix1;
			Minus(matrix3, matrix2);
		}
		static double Cofactor(std::vector<std::vector<double> > det, int R, int C) //���㷽��R��C�е�����ʽ
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
		static bool Adjoint(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2) //���㷽��İ�����
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
		static bool Inverse(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2) //���㷽�������
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
		static bool Multi(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> > matrix2, std::vector<std::vector<double> >& matrix3) //����������֮��
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
		static bool Power(std::vector<std::vector<double> > matrix1, std::vector<std::vector<double> >& matrix2, int power) //���㷽���power����
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
		static bool Eigen() //���㷽�������ֵ����������
		{

			return true;
		}

	protected:
		static void RowTrans(std::vector<std::vector<double> >& det, int r, double k) //��r�г�k��
		{
			int n = det[1].size() - 1;
			for (int i = 1; i <= n; i++) det[r][i] = k * det[r][i];
		}
		static void RowTrans(std::vector<std::vector<double> >& det, int r1, int r2, double k) //��r1�м���r2�е�k����r1 = r1 + k * r2
		{
			int n = det[1].size() - 1;
			for (int i = 1; i <= n; i++) det[r1][i] = det[r1][i] + k * det[r2][i];
		}
		static void RowSwap(std::vector<std::vector<double> >& det, int r1, int r2) //��������
		{
			det[r1].swap(det[r2]);
		}
	};

};