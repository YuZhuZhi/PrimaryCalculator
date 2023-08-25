#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>

class Linear
{
public:
	template <typename T>
	class Polynomial //����ʽ
	{
	public:
		Polynomial()
		{
			_poly.clear();
		}
		Polynomial(int num)
		{
			_poly.resize(num);
		}
		Polynomial(const Polynomial& poly)
		{
			_poly = poly._poly;
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
		void Input(const int n) //ָ�������Ķ���ʽ����
		{
			this->_poly.clear();
			double temp = 0;
			for (int i = 0; i <= n; i++) {
				std::cin >> temp;
				this->_poly.push_back(temp);
			}
		}
		void Print() const
		{
			std::cout << *(this->_poly.end() - 1) << "x^" << this->Highest();
			for (int i = this->Highest() - 1; i >= 0; i--) {
				if (this->_poly[i] > 0) std::cout << ' + ' << this->_poly[i] << "x^" << i;
				else if (this->_poly[i] < 0) std::cout << ' - ' << -this->_poly[i] << "x^" << i;
			}
		}
		int Terms() const
		{
			return (this->_poly.size());
		}
		int Highest() const
		{
			return (this->_poly.size() - 1);
		}
		Polynomial& operator=(const Polynomial& poly)
		{
			this->_poly = poly._poly;
		}
		Polynomial& operator-()
		{
			for (auto vit = this->_poly.begin(); vit != this->_poly.end(); vit++) *vit = -(*vit);
			return (*this);
		}
		Polynomial operator+(const Polynomial& poly) const //������ʽ���
		{
			if (this->Highest() >= poly.Highest()) {
				Polynomial result(*this);

			}
			else if (this->Highest() < poly.Highest()) {
				Polynomial result(poly);

			}

		}
		Polynomial operator-(const Polynomial& poly) const //������ʽ���
		{
			if (this->Highest() >= poly.Highest()) {
				Polynomial result(*this);

			}
			else if (this->Highest() < poly.Highest()) {
				Polynomial result(poly);

			}

		}
		Polynomial operator*(const Polynomial& poly) const //������ʽ���
		{
			Polynomial result(this->_poly.size() * poly._poly.size());


		}

	protected:
		std::vector<T> _poly;

	};

	template <typename T>
	class Matrix //����ʽ���������
	{
	public:
		Matrix()
		{
			_matrix.clear();
		}
		Matrix(const int n)
		{
			_matrix.clear();
			_matrix.resize(n + 1, std::vector<T>(n + 1));
		}
		Matrix(const int R, const int C)
		{
			_matrix.clear();
			_matrix.resize(R + 1, std::vector<T>(C + 1));
		}
		Matrix(const Matrix& matrix)
		{
			_matrix = matrix._matrix;
		}
		Matrix(const std::vector<std::vector<T> >& matrix)
		{
			_matrix = matrix;
		}
		void Input(const int n) //��������n�ķ�������
		{
			this->_matrix.clear();
			T temp = 0;
			std::vector<T> vectp;
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
		void Input(const int R, const int C) //��������R������C�ľ�������
		{
			this->_matrix.clear();
			T temp = 0;
			std::vector<T> vectp;
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
		void Print() const //��������
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
		int Rank() const //���ط���Ľ���
		{
			if (Row(*this) != Column(*this)) return 0;
			else return (this->_matrix.size() - 1);
		}
		int Row() const //���ؾ��������
		{
			return (this->_matrix.size() - 1);
		}
		int Column() const //���ؾ��������
		{
			return (this->_matrix[1].size() - 1);
		}
		bool SetRow(const std::vector<T>& target, const int r)
		{
			if (this->_matrix.size() == target.size()) {
				this->_matrix[r] = target;
				return true;
			}
			else return false;
		}
		bool SetColumn(const std::vector<T>& target, const int c)
		{
			if (this->_matrix[1].size() == target.size()) {
				this->Transpose();
				this->_matrix[c] = target;
				this->Transpose();
				return true;
			}
			else return false;
		}
		void SetElmn(const T target, const int r, const int c)
		{
			this->_matrix[r][c] = target;
		}
		T GetElmn(const int r, const int c)
		{
			return (this->_matrix[r][c]);
		}
		Matrix UnitMatrix(const int n) //���ɵ�λ����
		{
			Matrix matrix(n);
			for (int i = 1; i <= n; i++) matrix._matrix[i][i] = (T)1;
			return matrix;
		}
		Matrix UnitMatrix(const int n, const T multi) //���ɵ�λ�����multi��
		{
			Matrix matrix(n);
			for (int i = 1; i <= n; i++) this->_matrix[i][i] = multi;
			return matrix;
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
		Matrix& Transpose() //��ԭ����ת��
		{
			Transpose(*this);
			return (*this);
		}
		Matrix& AddRight(const Matrix& matrix) //��ԭ�����ұ߸���ͬ��������
		{
			if (AddRight(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix& AddDown(const Matrix& matrix) //��ԭ�������渽��ͬ��������
		{
			if (AddDown(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix& Plus(const Matrix& matrix) //��this����matirx
		{
			if (Plus(*this, matrix)) return (*this);
			else return (*this);
		}
		Matrix& Minus(const Matrix& matrix) //��this��ȥmatirx
		{
			if (Minus(*this, matrix)) return (*this);
			else return (*this);
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
		Matrix& Multi(const Matrix& matrix1, const Matrix& matrix2) //����matrix1��matrix2֮��
		{
			if (Multi(matrix1, matrix2, *this)) return (*this);
			else return (*this);
		}
		Matrix Power(const int power)
		{
			Matrix temp = *this;
			if (Power(*this, temp, power)) return (temp);
			else return (temp);
		}
		Matrix& operator=(const Matrix& matrix)
		{
			this->_matrix = matrix._matrix;
			return (*this);
		}
		Matrix operator+(const Matrix& matrix) const //���ı�this��matrix
		{
			int R = Row(*this), C = Column(*this);
			Matrix temp(R, C);
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp._matrix[r][c] = this->_matrix[r][c] + matrix._matrix[r][c];
			}
			return temp;
		}
		Matrix operator-(const Matrix& matrix) const //���ı�this��matrix
		{
			int R = Row(*this), C = Column(*this);
			Matrix temp(R, C);
			for (int r = 1; r <= R; r++) {
				for (int c = 1; c <= C; c++) temp._matrix[r][c] = this->_matrix[r][c] - matrix._matrix[r][c];
			}
			return temp;
		}
		Matrix operator*(const Matrix& matrix)  //���ı�this��matrix
		{
			Matrix result;
			if (Multi(*this, matrix, result)) return result;
			else return result;
		}
		Matrix operator/(const Matrix& matrix) 
		{
			Matrix temp = matrix, result;
			if (Multi(*this, temp.Inverse(), result)) return result;
			else return result;
		}
		bool Eigen() //���㷽�������ֵ����������
		{

			return true;
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
		static T Det(const Matrix& det) //��������ʽ��ֵ
		{
			Matrix temp = det;
			int n = Rank(temp);
			if (n == 0) return -1;
			temp = temp.GausElmn(); //��˹��ԪΪ�Ͻ�����
			T multi = 1;
			for (int i = 1; i <= n; i++) multi = multi * temp._matrix[i][i];
			return (multi);
		}
		static T Cofactor(const Matrix& det, const int R, const int C) //���㷽��R��C�е�����ʽ
		{
			Matrix temp(Row(det) - 1, Column(det) - 1);
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

	protected:
		std::vector<std::vector<T> > _matrix;
		void RowTrans(Matrix& matrix, const int r, const T k) //��r�г�k��
		{
			int n = matrix._matrix[1].size() - 1;
			for (int i = 1; i <= n; i++) matrix._matrix[r][i] = k * matrix._matrix[r][i];
		}
		void RowTrans(Matrix& matrix, const int r1, const int r2, const T k) //��r1�м���r2�е�k����r1 = r1 + k * r2
		{
			int n = matrix._matrix[1].size() - 1;
			for (int i = 1; i <= n; i++) matrix._matrix[r1][i] = matrix._matrix[r1][i] + k * matrix._matrix[r2][i];
		}
		void RowSwap(Matrix& matrix, const int r1, const int r2) //��������
		{
			matrix._matrix[r1].swap(matrix._matrix[r2]);
		}
		void GausElmn(Matrix& matrix) //��ԭ�����˹��Ԫ�������н�����
		{
			int R = Row(matrix), C = Column(matrix);
			for (int c = 1; c <= C; c++) {
				int r = 1;
				for (r = c; r <= R; r++) {
					if (fabs((double)matrix._matrix[r][c]) >= 1e-6) break;
				}
				if (r <= R) {
					int temp = r;
					RowSwap(matrix, r, c);
					for (r = temp + 1; r <= R; r++) {
						if (fabs((double)matrix._matrix[r][c]) <= 1e-6) break;
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
			Matrix temp(R, C);
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
			T M = Det(matrix1);
			if (Row(matrix1) != Column(matrix1) || M == (T)0) return false;
			else {
				for (int r = 1; r <= Rank(matrix1); r++) {
					for (int c = 1; c <= Rank(matrix1); c++) matrix2._matrix[r][c] = ((T)Cofactor(matrix1, c, r) / M) * (T)pow(-1, r + c);
				}
				return true;
			}
		}
		bool Multi(const Matrix& matrix1, const Matrix& matrix2, Matrix& matrix3) //����������֮�����������matrix3
		{
			int R = Row(matrix1), C = Column(matrix2);
			matrix3._matrix.clear();
			matrix3._matrix.resize(matrix1._matrix.size(), std::vector<T>(matrix2._matrix[1].size(), 0));
			if (matrix1._matrix[1].size() != matrix2._matrix.size()) return false;
			else {
				Matrix inv_matrix2 = matrix2;
				Transpose(inv_matrix2);
				for (int r = 1; r <= R; r++) {
					for (int c = 1; c <= C; c++) {
						matrix3._matrix[r][c] = VectorMulti(matrix1._matrix[r], inv_matrix2._matrix[c]);
					}
				}
				return true;
			}
		}
		bool Power(const Matrix& matrix1, Matrix& matrix2, const int power) //���㷽���power����
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
		T VectorMulti(const std::vector<T>& A, const std::vector<T>& B)
		{
			T result = 0;
			for (int i = 1; i <= A.size() - 1; i++) result = result + A[i] * B[i];
			return result;
		}

	};

	template <typename T>
	class Vector //��������
	{
	public:
		Vector()
		{
			_vector.clear();
		}
		Vector(const int n)
		{
			_vector.clear();
			_vector.resize(n + 1);
		}
		Vector(const Vector& vector)
		{
			_vector = vector._vector;
		}
		Vector(const std::vector<double>& vector)
		{
			_vector = vector;
		}
		void Input() //��������
		{
			this->_vector.clear();
			this->_vector.push_back(0); //0λ����0���
			T temp = 0;
			while (std::cin >> temp) { //��δ����Ctrl+Zǰ�����ǽ�������
				this->_vector.push_back(temp);
			}
			std::cin.clear(); //���Ctrl+Z״̬
		}
		void Input(const int n) //ָ��ά������������
		{
			this->_vector.clear();
			this->_vector.push_back(0); //0λ����0���
			T temp = 0;
			for (int i = 1; i <= n; i++) {
				std::cin >> temp;
				this->_vector.push_back(temp);
			}
		}
		void Print() const //�������
		{
			auto temp = this->_vector;
			std::cout << "(";
			for (auto vdit = temp.begin() + 1; vdit != temp.end() - 1; vdit++) {
				std::cout << *vdit << ", ";
			}
			std::cout << *(this->_vector.end() - 1) << ")" << std::endl << std::endl;
		}
		int Dimension() const
		{
			return (this->_vector.size() - 1);
		}
		Vector& Plus(const Vector& vector)
		{
			for (int i = 1; i <= this->Dimension(); i++) this->_vector[i] = this->_vector[i] + vector._vector[i];
			return (*this);
		}
		Vector& Minus(const Vector& vector)
		{
			for (int i = 1; i <= this->Dimension(); i++) this->_vector[i] = this->_vector[i] - vector._vector[i];
			return (*this);
		}
		Vector VectorPro(const Vector& vector) 
		{
			Vector temp;
			if (VectorPro(*this, vector, temp)) return (temp);
			else return (temp);
		}
		double ScalarPro(const Vector& vector) const
		{
			return (ScalarPro(*this, vector));
		}
		double Length() const //���������ĳ���
		{
			return (sqrt(ScalarPro(*this, *this)));
		}
		double Norm(const int power) const
		{
			int n = this->Dimension();
			double temp = 0;
			for (int i = 1; i <= n; i++) {
				temp = temp + pow(this->_vector[i], power);
			}
			return (pow(temp, 1.0 / (double)power));
		}
		double Angle(const Vector& vector) const
		{
			return Angle(*this, vector);
		}

		static void Print(const Vector& vector)
		{
			Vector temp = vector;
			std::cout << "(";
			for (auto vdit = temp._vector.begin() + 1; vdit != temp._vector.end() - 1; vdit++) {
				std::cout << *vdit << ", ";
			}
			std::cout << *(temp._vector.end() - 1) << ")" << std::endl << std::endl;
		}
		static double Dimension(const Vector& vector)
		{
			return (vector._vector.size() - 1);
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
		static double Length(const Vector x) //���������ĳ���
		{
			return (sqrt(ScalarPro(x, x))); //���������ݼ���ƽ����
		}
		static double Angle(const Vector& x, const Vector& y) //������������н�
		{
			return (acos((ScalarPro(x, y)) / (Length(x) * Length(y))));
		}
		static double Norm(const Vector& x, const int power)
		{
			int n = Dimension(x);
			double temp = 0;
			for (int i = 1; i <= n; i++) {
				temp = temp + pow(x._vector[i], power);
			}
			return (pow(temp, 1.0 / power));
		}

	protected:
		std::vector<T> _vector;
		bool VectorPro(const Vector& x, const Vector& y, Vector& result) //������������������result������������
		{
			if (Dimension(x) == 3 && Dimension(y) == 3) {
				result._vector.clear();
				result._vector.push_back(0);
				Matrix<T> temp(3, 3); //����һ����������ʽ
				std::vector<double> xtemp = x._vector, ytemp = y._vector;
				temp.SetRow(xtemp, 2);
				temp.SetRow(ytemp, 3);
				for (int i = 1; i <= 3; i++) result._vector.push_back(Matrix<T>::Cofactor(temp, 1, i)); //�����Ӧ����ʽ
				result._vector[2] = -result._vector[2]; //�ڶ�ά��ֵȡ�෴��
				return true;
			}
			else return false;
		}

	};

	class Equation //�ⷽ��
	{
	public:
		static void Print(std::vector<double>& result) //������̵Ķ����
		{
			int n = result.size() - 1;
			for (int i = 1; i <= n; i++) {
				std::cout << "x_" << i << " = " << result[i] << std::endl;
			}
		}
		static bool Quadratic(const double a, const double b, const double c, std::vector<double>& result) //����η��̣���result��װ��
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
		static bool Cubic(const double a, const double b, const double c, const double d, std::vector<double>& result) //�����η��̣���result��װ��
		{
			return true;
		}
		static bool Set(const Matrix<double> coeff, const std::vector<double> right, std::vector<double>& result) //�ⷽ���飬��result��װ��
		{
			result.clear();
			result.push_back(Matrix<double>::Det(coeff)); //0λ������ϵ������ʽ
			if (result[0] == 0) return false;
			else {
				for (int i = 1; i <= Matrix<double>::Rank(coeff); i++) {
					Matrix<double> temp = coeff;
					temp.Transpose();
					temp.SetRow(right, i);
					result.push_back((Matrix<double>::Det(temp)) / result[0]);
				}
				return true;
			}
		}
	};

};