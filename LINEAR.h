#pragma once
#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#define Max(a, b) (a > b) ? a : b;
#define Min(a, b) (a > b) ? b : a;

class Linear
{
public:
	class Polynomial //����ʽ
	{
	public:
		Polynomial()
		{
			_poly.clear();
		}
		void Input() //����ʽϵ������
		{
			this->_poly.clear();
			double temp = 0;
			while (std::cin >> temp) { //��δ����Ctrl+Zǰ�����ǽ�������
				this->_poly.push_back(temp);
			}
			std::cin.clear(); //���Ctrl+Z״̬
		}
		void Input(int n) //ָ�������Ķ���ʽ����
		{
			this->_poly.clear();
			double temp = 0;
			for (int i = 0; i <= n; i++) {
				std::cin >> temp;
				this->_poly.push_back(temp);
			}
		}

		static Polynomial Multi(const Polynomial& A, const Polynomial& B, Polynomial& C) //������ʽ���
		{
			int max = Max(A._poly.size(), B._poly.size());
			C._poly.resize(max);

		}

	protected:
		std::vector<double> _poly;

	};

	class Matrix //����ʽ���������
	{
	public:
		Matrix()
		{
			_matrix.clear();
		}
		Matrix(int n)
		{
			_matrix.clear();
			_matrix.resize(n + 1, std::vector<double>(n + 1));
		}
		Matrix(int R, int C)
		{
			_matrix.clear();
			_matrix.resize(R + 1, std::vector<double>(C + 1));
		}
		void Input(int n) //��������n�ķ�������
		{
			this->_matrix.clear();
			double temp = 0;
			std::vector<double> vectp;
			for (int i = 0; i <= n; i++) vectp.push_back(0);
			this->_matrix.push_back(vectp); //��άvector��0������0Ԫ��
			for (int r = 1; r <= n; r++) {
				vectp.resize(1, 0);
				for (int c = 1; c <= n; c++) {
					std::cin >> temp;
					vectp.push_back(temp);
				}
				this->_matrix.push_back(vectp);
			}
		}
		void Input(int R, int C) //��������R������C�ľ�������
		{
			this->_matrix.clear();
			double temp = 0;
			std::vector<double> vectp;
			for (int i = 1; i <= C + 1; i++) vectp.push_back(0);
			this->_matrix.push_back(vectp); //��άvector��0������0Ԫ��
			for (int r = 1; r <= R; r++) {
				vectp.resize(1, 0);
				for (int c = 1; c <= C; c++) {
					std::cin >> temp;
					vectp.push_back(temp);
				}
				this->_matrix.push_back(vectp);
			}
		}
		void Print() //��������
		{
			for (int r = 1; r <= Row(*this); r++) {
				for (int c = 1; c <= Column(*this); c++) {
					std::cout.width(6); //�ܼ���Ԫ�ؽ�С��ÿ��Ԫ����6��Ŀռ����
					std::cout.flags(std::ios::right); //�̶������λС�������Ҷ���
					std::cout.precision(2); //�̶������λС��
					std::cout << std::fixed << this->_matrix[r][c] << " ";
				}
				std::cout << std::endl << std::endl;
			}
		}
		static void Print(const Matrix& matrix) //��������
		{
			for (int r = 1; r <= Row(matrix); r++) {
				for (int c = 1; c <= Column(matrix); c++) {
					std::cout.width(6); //�ܼ���Ԫ�ؽ�С��ÿ��Ԫ����6��Ŀռ����
					std::cout.flags(std::ios::right); //�̶������λС�������Ҷ���
					std::cout.precision(2); //�̶������λС��
					std::cout << std::fixed << matrix._matrix[r][c] << " ";
				}
				std::cout << std::endl << std::endl;
			}
		}
		int Rank() //���ط���Ľ���
		{
			if (Row(*this) != Column(*this)) return 0;
			else return (this->_matrix.size() - 1);
		}
		int Row() //���ؾ��������
		{
			return (Row(*this));
		}
		int Column() //���ؾ��������
		{
			return (Column(*this));
		}
		bool ChangeRow(std::vector<double> target, int r)
		{
			if (this->_matrix.size() == target.size()) {
				this->_matrix[r] = target;
				return true;
			}
			else return false;
		}
		bool ChangeColumn(std::vector<double> target, int c)
		{
			if (this->_matrix[1].size() == target.size()) {
				this->Transpose();
				this->_matrix[c] = target;
				this->Transpose();
				return true;
			}
			else return false;
		}
		void ChangeElmn(double target, int r, int c)
		{
			this->_matrix[r][c] = target;
		}
		Matrix UnitMatrix(int n) //���ɵ�λ����
		{
			this->_matrix.clear();
			this->_matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) this->_matrix[i][i] = 1;
			return (*this);
		}
		Matrix UnitMatrix(int n, int num) //���ɵ�λ�����num��
		{
			this->_matrix.clear();
			this->_matrix.resize(n + 1, std::vector<double>(n + 1, 0));
			for (int i = 1; i <= n; i++) this->_matrix[i][i] = num;
			return (*this);
		}
		Matrix GausElmn() //����ԭ�����˹��Ԫ������н�����
		{
			Matrix temp = *this;
			GausElmn(temp);
			return (temp);
		}
		Matrix RowSim() //����ԭ������������
		{
			Matrix temp = *this;
			RowSim(temp);
			return (temp);
		}
		static double Det(const Matrix& det) //��������ʽ��ֵ
		{
			Matrix temp = det;
			int n = Rank(temp);
			if (n == 0) return -1;
			temp = temp.GausElmn(); //��˹��ԪΪ�Ͻ�����
			double multi = 1;
			for (int i = 1; i <= n; i++) multi = multi * temp._matrix[i][i];
			return (multi);
		}
		Matrix Transpose() //��ԭ����ת��
		{
			Transpose(*this);
			return (*this);
		}
		Matrix AddRight(const Matrix& matrix) //��ԭ�����ұ߸���ͬ��������
		{
			if (AddRight(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix AddDown(const Matrix& matrix) //��ԭ�������渽��ͬ��������
		{
			if (AddDown(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix Plus(const Matrix& matrix) //��this����matirx
		{
			if (Plus(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix Minus(const Matrix& matrix) //��this��ȥmatirx
		{
			if (Minus(*this, matrix)) return (*this);
			else return (*this);
		}
		static double Cofactor(const Matrix& det, int R, int C) //���㷽��R��C�е�����ʽ
		{
			Matrix temp(R - 1, C - 1);
			for (int r = 1; r <= Row(det); r++) {
				if (r < R) {
					for (int c = 1; c <= Column(det); c++) {
						if (c < C) temp._matrix[r][c] = det._matrix[r][c];
						if (c > C) temp._matrix[r][c - 1] = det._matrix[r][c];
					}
				}
				if (r > R) {
					for (int c = 1; c <= Column(det); c++) {
						if (c < C) temp._matrix[r - 1][c] = det._matrix[r][c];
						if (c > C) temp._matrix[r - 1][c - 1] = det._matrix[r][c];
					}
				}
			}
			return (Det(temp));
		}
		Matrix Adjoint() //���ط���matrix�İ�����
		{
			Matrix temp = *this;
			int n = Rank(temp);
			for (int r = 1; r <= n; r++) {
				for (int c = 1; c <= n; c++) temp._matrix[r][c] = Cofactor(*this, c, r) * pow(-1, r + c);
			}
			return (temp);
		}
		Matrix Inverse() //����ԭ���������
		{
			Matrix temp = *this;
			if (Inverse(*this, temp)) return temp;
			else return temp;
		}
		Matrix Multi(const Matrix& matrix1, const Matrix& matrix2) //����matrix1��matrix2֮��
		{
			if (Multi(matrix1, matrix2, *this)) {
				return (*this);
			}
			else {
				return (*this);
			}
		}
		Matrix Power(int power)
		{
			Matrix temp = *this;
			if (Power(*this, temp, power)) return (*this);
			else return (*this);
		}
		bool Eigen() //���㷽�������ֵ����������
		{

			return true;
		}
		Matrix operator=(const Matrix& matrix)
		{
			this->_matrix = matrix._matrix;
			return (*this);
		}
		Matrix operator+(const Matrix& matrix) //���ı�this��matrix
		{
			int R = Row(*this), C = Column(*this);
			Matrix temp(R, C);
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp._matrix[r][c] = this->_matrix[r][c] + matrix._matrix[r][c];
			}
			return temp;
		}
		Matrix operator-(const Matrix& matrix) //���ı�this��matrix
		{
			int R = Row(*this), C = Column(*this);
			Matrix temp(R, C);
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp._matrix[r][c] = this->_matrix[r][c] - matrix._matrix[r][c];
			}
			return temp;
		}
		Matrix operator*(const Matrix& matrix) //���ı�this��matrix
		{
			Matrix temp;
			if (Multi(*this, matrix, temp)) return temp;
			else return temp;
		}


	protected:
		std::vector<std::vector<double> > _matrix;
		static int Rank(const Matrix& matrix) //���ط���Ľ���
		{
			if (Row(matrix) != Column(matrix)) return 0;
			else return (matrix._matrix.size() - 1);
		}
		static int Row(const Matrix& matrix) //���ؾ��������
		{
			return (matrix._matrix.size() - 1);
		}
		static int Column(const Matrix& matrix) //���ؾ��������
		{
			return (matrix._matrix[1].size() - 1);
		}
		void RowTrans(Matrix& matrix, int r, double k) //��r�г�k��
		{
			int n = matrix._matrix[1].size() - 1;
			for (int i = 1; i <= n; i++) matrix._matrix[r][i] = k * matrix._matrix[r][i];
		}
		void RowTrans(Matrix& matrix, int r1, int r2, double k) //��r1�м���r2�е�k����r1 = r1 + k * r2
		{
			int n = matrix._matrix[1].size() - 1;
			for (int i = 1; i <= n; i++) matrix._matrix[r1][i] = matrix._matrix[r1][i] + k * matrix._matrix[r2][i];
		}
		void RowSwap(Matrix& matrix, int r1, int r2) //��������
		{
			matrix._matrix[r1].swap(matrix._matrix[r2]);
		}
		void GausElmn(Matrix& matrix) //��ԭ�����˹��Ԫ�������н�����
		{
			int R = Row(matrix), C = Column(matrix);
			for (int c = 1; c <= C; c++) {
				int r = 1;
				for (r = c; r <= R; r++) {
					if (fabs(matrix._matrix[r][c]) >= 1e-6) break;
				}
				if (r <= R) {
					int temp = r;
					RowSwap(matrix, r, c);
					for (r = temp + 1; r <= R; r++) {
						if (matrix._matrix[r][c] == 0) {}
						else RowTrans(matrix, r, temp, -(matrix._matrix[r][c] / matrix._matrix[temp][c]));
					}
				}
			}
		}
		void RowSim(Matrix& matrix) //��ԭ���������������
		{
			int R = Row(matrix), C = Column(matrix);
			GausElmn(matrix); //��Ϊ�н�����
			for (int r = 1; r <= R; r++) {
				int c = 1;
				for (c = 1; c <= C; c++) {
					if (fabs(matrix._matrix[r][c]) >= 1e-6) break;
				}
				RowTrans(matrix, r, 1.0 / matrix._matrix[r][c]);
			}
			for (int r = R; r >= 1; r--) {
				int c = C;
				for (c = 1; c <= C; c++) {
					if (fabs(matrix._matrix[r][c]) >= 1e-6) break;
				}
				for (int r2 = r - 1; r2 >= 1; r2--) {
					RowTrans(matrix, r2, r, -matrix._matrix[r2][c]);
				}
			}
		}
		void Transpose(Matrix& matrix) //��ԭ����ת��
		{
			int R = Column(matrix), C = Row(matrix);
			Matrix temp(C, R);
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp._matrix[r][c] = matrix._matrix[c][r];
			}
			matrix = temp;
		}
		bool AddRight(Matrix& matrix1, const Matrix& matrix2) //��ԭ�����ұ߸���ͬ��������
		{
			if (Row(matrix1) != Row(matrix2)) return false;
			else {
				for (int r = 1; r <= Row(matrix1); r++) {
					std::vector<double> temp = matrix1._matrix[r];
					temp.insert(temp.end(), matrix2._matrix[r].begin() + 1, matrix2._matrix[r].end());
					matrix1._matrix[r] = temp;
				}
				return true;
			}
		}
		bool AddDown(Matrix& matrix1, const Matrix& matrix2) //��ԭ�������渽��ͬ��������
		{
			if (matrix1._matrix[1].size() != matrix2._matrix[1].size()) return false;
			else {
				for (int r = 1; r <= Row(matrix2); r++) matrix1._matrix.push_back(matrix2._matrix[r]);
				return true;
			}
		}
		bool Plus(Matrix& matrix1, const Matrix& matrix2) //��matrix1����matirx2
		{
			if (Row(matrix1) == Row(matrix2) && Column(matrix1) == Column(matrix2)) {
				int R = Row(matrix1), C = Column(matrix1);
				for (int r = 1; r <= R; r++) {
					for (int c = 1; c <= C; c++) matrix1._matrix[r][c] = matrix1._matrix[r][c] + matrix2._matrix[r][c];
				}
				return true;
			}
			else return false;
		}
		bool Minus(Matrix& matrix1, const Matrix& matrix2) //��matrix1��ȥmatirx2
		{
			if (Row(matrix1) == Row(matrix2) && Column(matrix1) == Column(matrix2)) {
				int R = Row(matrix1), C = Column(matrix1);
				for (int r = 1; r <= R; r++) {
					for (int c = 1; c <= C; c++) matrix1._matrix[r][c] = matrix1._matrix[r][c] - matrix2._matrix[r][c];
				}
				return true;
			}
			else return false;
		}
		bool Inverse(const Matrix& matrix1, Matrix& matrix2) //���㷽������󣬽������matrix2
		{
			matrix2 = matrix1;
			double M = Det(matrix1);
			if (Row(matrix1) != Column(matrix1) || M == 0) return false;
			else {
				for (int r = 1; r <= Rank(matrix1); r++) {
					for (int c = 1; c <= Rank(matrix1); c++) matrix2._matrix[r][c] = (Cofactor(matrix1, c, r) / M) * pow(-1, r + c);
				}
				return true;
			}
		}
		bool Multi(const Matrix& matrix1, const Matrix& matrix2, Matrix& matrix3) //����������֮�����������matrix3
		{
			int R = Row(matrix1), C = Column(matrix2);
			matrix3._matrix.clear();
			matrix3._matrix.resize(matrix1._matrix.size(), std::vector<double>(matrix2._matrix[1].size(), 0));
			if (matrix1._matrix[1].size() != matrix2._matrix.size()) return false;
			else {
				Matrix inv_matrix2 = matrix2;
				Transpose(inv_matrix2);
				for (int r = 1; r <= R; r++) {
					for (int c = 1; c <= C; c++) {
						Vector x(matrix1._matrix[r]), y(inv_matrix2._matrix[c]);
						matrix3._matrix[r][c] = Vector::ScalarPro(x, y);
					}
				}
				return true;
			}
		}
		bool Power(const Matrix& matrix1, Matrix& matrix2, int power) //���㷽���power����
		{
			if (Row(matrix1) != Column(matrix1)) return false;
			else {
				matrix2 = matrix1;
				if (power > 1) {
					for (int i = 2; i <= power; i++) {
						matrix2 = matrix2 * matrix1;
					}
				}
				if (power == 1) return true;
				if (power == 0) {
					Matrix temp = UnitMatrix(Row(matrix1));
					matrix2 = temp;
				}
				if (power < 0) {
					Matrix inv_matrix2;
					inv_matrix2 = matrix2.Inverse();
					matrix2 = inv_matrix2;
					for (int i = 2; i <= -power; i++) {
						matrix2 = matrix2 * inv_matrix2;
					}
				}
				return true;
			}
		}

	};

	class Vector //��������
	{
	public:
		Vector()
		{
			_vector.clear();
		}
		Vector(int n)
		{
			_vector.clear();
			_vector.resize(n + 1);
		}
		Vector(std::vector<double> vector)
		{
			_vector = vector;
		}
		void Input() //��������
		{
			this->_vector.clear();
			this->_vector.push_back(0); //0λ����0���
			double temp = 0;
			while (std::cin >> temp) { //��δ����Ctrl+Zǰ�����ǽ�������
				this->_vector.push_back(temp);
			}
			std::cin.clear(); //���Ctrl+Z״̬
		}
		void Input(int n) //ָ��ά������������
		{
			this->_vector.clear();
			this->_vector.push_back(0); //0λ����0���
			double temp = 0;
			for (int i = 1; i <= n; i++) {
				std::cin >> temp;
				this->_vector.push_back(temp);
			}
		}
		void Print() //�������
		{
			std::cout << "(";
			for (std::vector<double>::iterator vdit = this->_vector.begin() + 1; vdit != this->_vector.end() - 1; vdit++) {
				std::cout << *vdit << ", ";
			}
			std::cout << *(this->_vector.end() - 1) << ")" << std::endl << std::endl;
		}
		static void Print(Vector& vector)
		{
			std::cout << "(";
			for (std::vector<double>::iterator vdit = vector._vector.begin() + 1; vdit != vector._vector.end() - 1; vdit++) {
				std::cout << *vdit << ", ";
			}
			std::cout << *(vector._vector.end() - 1) << ")" << std::endl << std::endl;
		}
		double Dimension()
		{
			return (this->_vector.size() - 1);
		}
		Vector Plus(const Vector& vector)
		{
			for (int i = 1; i <= this->Dimension(); i++) this->_vector[i] = this->_vector[i] + vector._vector[i];
			return (*this);
		}
		Vector Minus(const Vector& vector)
		{
			for (int i = 1; i <= this->Dimension(); i++) this->_vector[i] = this->_vector[i] - vector._vector[i];
			return (*this);
		}
		static double ScalarPro(const Vector& x, const Vector& y) //����������������
		{
			if (Dimension(x) != Dimension(y)) return (0);
			double crssum = 0;
			for (int i = 1; i <= Dimension(x); i++) {
				crssum = crssum + x._vector[i] * y._vector[i];
			}
			return (crssum);
		}
		double ScalarPro(const Vector& vector)
		{
			return (ScalarPro(*this, vector));
		}
		Vector VectorPro(const Vector& vector)
		{
			Vector temp;
			if (VectorPro(*this, vector, temp)) return (temp);
			else return (temp);
		}
		double Length() //���������ĳ���
		{
			return (sqrt(ScalarPro(*this, *this))); //���������ݼ���ƽ����
		}
		static double Length(Vector x) //���������ĳ���
		{
			return (sqrt(ScalarPro(x, x))); //���������ݼ���ƽ����
		}
		static double Angle(const Vector& x, const Vector& y) //������������н�
		{
			return (acos((ScalarPro(x, y)) / (Length(x) * Length(y))));
		}

	protected:
		std::vector<double> _vector;
		static double Dimension(const Vector& vector)
		{
			return (vector._vector.size() - 1);
		}
		bool VectorPro(const Vector& x, const Vector& y, Vector& result) //������������������result������������
		{
			if (Dimension(x) == 3 && Dimension(y) == 3) {
				result._vector.clear();
				result._vector.push_back(0);
				Matrix temp(3, 3); //����һ����������ʽ
				std::vector<double> xtemp = x._vector, ytemp = y._vector;
				temp.ChangeRow(xtemp, 2);
				temp.ChangeRow(ytemp, 3);
				for (int i = 1; i <= 3; i++) result._vector.push_back(Matrix::Cofactor(temp, 1, i)); //�����Ӧ����ʽ
				result._vector[2] = -result._vector[2]; //�ڶ�ά��ֵȡ�෴��
				return true;
			}
			else return false;
		}

	};

	class Equation //�ⷽ��
	{
	public:
		static void Print(std::vector<double> result) //������̵Ķ����
		{
			int n = result.size() - 1;
			for (int i = 1; i <= n; i++) {
				std::cout << "x_" << i << " = " << result[i] << std::endl;
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
		static bool Set(Matrix coeff, std::vector<double> right, std::vector<double>& result) //�ⷽ���飬��result��װ��
		{
			result.clear();
			result.push_back(Matrix::Det(coeff)); //0λ������ϵ������ʽ
			if (result[0] == 0) return false;
			else {
				for (int i = 1; i <= coeff.Rank(); i++) {
					Matrix temp = coeff;
					temp.Transpose();
					temp.ChangeRow(right, i);
					result.push_back((Matrix::Det(temp)) / result[0]);
				}
				return true;
			}
		}
	};

//protected:
//	static double ScalarPro(std::vector<double> x, std::vector<double> y)
//	{
//		double 
//	}

};