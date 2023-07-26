#pragma once
#include <iostream>
#include "DATA.h"

extern char YN;
extern int choicefir, choicesec, choicethr;
extern std::vector<double> data_set, x_set, y_set, right_value;
extern Linear::Vector x, y;
extern Linear::Matrix matrix, matrix1, matrix2, matrix3, coeff;

class Operate
{
public:
	class CF1
	{
	public:
		static void CS1()
		{
			if (YN == 'N' || YN == 'n') {
				std::cout << "请输入数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
				Data::Process::Input(data_set);
				std::cout << std::endl << "总和是：" << std::endl << Data::Process::Summary(data_set) << std::endl;
			}
			if (YN == 'Y' || YN == 'y') {
				std::cout << std::endl << "总和是：" << std::endl << Data::Process::Summary(data_set) << std::endl;
			}
		}
		static void CS2()
		{
			if (YN == 'N' || YN == 'n') {
				std::cout << "请输入数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
				Data::Process::Input(data_set);
				std::cout << std::endl << "平均数是：" << std::endl << Data::Process::Average(data_set) << std::endl;
			}
			if (YN == 'Y' || YN == 'y') {
				std::cout << std::endl << "平均数是：" << std::endl << Data::Process::Average(data_set) << std::endl;
			}
		}
		static void CS3()
		{
			if (YN == 'N' || YN == 'n') {
				std::cout << "请输入数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
				Data::Process::Input(data_set);
				std::cout << std::endl << "标准差(数) / 实验标准差(物)是：" << std::endl << Data::Process::StdDev(data_set) << std::endl;
			}
			if (YN == 'Y' || YN == 'y') {
				std::cout << std::endl << "标准差(数) / 实验标准差(物)是：" << std::endl << Data::Process::StdDev(data_set) << std::endl;
			}
		}
		static void CS4()
		{
			if (YN == 'N' || YN == 'n') {
				std::cout << "请输入数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
				Data::Process::Input(data_set);
				std::cout << std::endl << "平均值标准差(数) / 平均值的实验标准差(物)是：" << std::endl << Data::Process::AverStdDev(data_set) << std::endl;
			}
			if (YN == 'Y' || YN == 'y') {
				std::cout << std::endl << "平均值标准差(数) / 平均值的实验标准差(物)是：" << std::endl << Data::Process::AverStdDev(data_set) << std::endl;
			}
		}
		static void CS5()
		{
			double S_A = 0, S_B = 0, S = 0, vamin = 0, v_A = 0, v_B = 1, v = 0, t_p = 0;
			int type = 0;
			char free = 'Y';
			if (YN == 'N' || YN == 'n') {
				std::cout << "请输入数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
				Data::Process::Input(data_set);
				if (Data::Process::Capacity(data_set) > 1) {
					S_A = Data::Process::AverStdDev(data_set);
					v_A = Data::Process::Capacity(data_set) - 1;
				}
				else {
					S_A = 0;
					v_A = 1;
				}
				std::cout << "请输入测量仪器的最小分度值：" << std::endl;
				std::cin >> vamin;
				std::cout << "请输入误差分布类型：" << std::endl;
				std::cout << "1\t均匀分布" << std::endl;
				std::cout << "2\t三角分布" << std::endl;
				std::cout << "3\t正态分布" << std::endl;
				std::cin >> type;
				switch (type) {
				case 1: 
					S_B = vamin / (sqrt(3));
					break;
				case 2: 
					S_B = vamin / 2.0;
					break;
				case 3: 
					S_B = vamin / (sqrt(6));
					break;
				default:
					std::cout << "ERROR!" << std::endl;
					CS5();
					break;
				}
				S = sqrt(S_A * S_A + S_B * S_B);
				std::cout << "S_A = " << S_A << std::endl;
				std::cout << "S_B = " << S_B << std::endl;
				std::cout << "S = " << S << std::endl << std::endl;
				std::cout << "自由度是否取无限大？(置信概率默认取0.95)" << std::endl;
				std::cout << "Y\t是" << std::endl;
				std::cout << "N\t否" << std::endl;
				std::cin >> free;
				if (free == 'Y' || free == 'y') {
					std::cout << "拓展不确定度delta N为：" << std::endl << S * 1.96 << std::endl;
				}
				if (free == 'N' || free == 'n') {
					double t_p = 1.96;
					v = pow(S, 4.0) / ((pow(S_A, 4.0) / v_A) + (pow(S_B, 4.0) / v_B));
					std::cout << "有效自由度v为：" << std::endl << v << std::endl;
					std::cout << "查询后输入置信系数t_p：" << std::endl;
					std::cin >> t_p;
					std::cout << "拓展不确定度delta N为：" << std::endl << S * t_p << std::endl;
				}
			}
			if (YN == 'Y' || YN == 'y') {
				std::cout << "S_A = " << S_A << std::endl;
				std::cout << "S_B = " << S_B << std::endl;
				std::cout << "S = " << S << std::endl << std::endl;
				std::cout << "自由度是否取无限大？(置信概率默认取0.95)" << std::endl;
				std::cout << "Y\t是" << std::endl;
				std::cout << "N\t否" << std::endl;
				std::cin >> free;
				if (free == 'Y' || free == 'y') {
					std::cout << "拓展不确定度delta N为：" << std::endl << S * 1.96 << std::endl;
				}
				if (free == 'N' || free == 'n') {
					double t_p = 1.96;
					v = pow(S, 4.0) / ((pow(S_A, 4.0) / v_A) + (pow(S_B, 4.0) / v_B));
					std::cout << "有效自由度v为：" << std::endl << v << std::endl;
					std::cout << "查询后输入置信系数t_p：" << std::endl;
					std::cin >> t_p;
					std::cout << "拓展不确定度delta N为：" << std::endl << S * t_p << std::endl;
				}
			}
		}
		static void CS6()
		{
			std::cout << "拟合方程为y = Kx + B" << std::endl;
			std::cout << "请输入自变量数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
			Data::Process::Input(x_set);
			std::cout << "请输入因变量数据：(键入Ctrl+Z后Enter以确认结束输入)" << std::endl;
			Data::Process::Input(y_set);
			std::cout << "x平均值是:" << std::endl << Data::Process::Average(x_set) << std::endl;
			std::cout << "y平均值是:" << std::endl << Data::Process::Average(y_set) << std::endl;
			std::cout << "x^2平均值是:" << std::endl << Data::Process::SquareAver(x_set) << std::endl;
			std::cout << "y^2平均值是:" << std::endl << Data::Process::SquareAver(y_set) << std::endl;
			std::cout << "xy平均值是:" << std::endl << Data::Process::CroseAver(x_set, y_set) << std::endl;
			std::cout << "拟合结果为：" << std::endl << "K = " << Data::LinearFit::K(x_set, y_set) << std::endl << "B = " << Data::LinearFit::B(x_set, y_set) << std::endl;
			std::cout << "相关系数R = " << Data::LinearFit::R(x_set, y_set) << std::endl;
		}
	};

	class CF2
	{
	public:
		static void CS1()
		{

		}
		static void CS2()
		{

		}
		static void CS3()
		{

		}
		static void CS4()
		{

		}
		static void CS5()
		{

		}
	};

	class CF3
	{
	public:
		class CS1
		{
		public:
			static void CT1()
			{
				double a = 0, b = 0, c = 0;
				std::vector<double> result;
				std::cout << "方程形式为ax^2 + bx + c = 0 (a不为0)" << std::endl;
				std::cout << "请输入系数a, b, c：(无需输入Ctrl+Z以确认)" << std::endl;
				std::cin >> a >> b >> c;
				Linear::Equation::Quadratic(a, b, c, result);
				Linear::Equation::Print(result);
			}
			static void CT2()
			{
				double a = 0, b = 0, c = 0, d = 0;
				std::vector<double> result;
				std::cout << "方程形式为ax^3 + bx^2 + cx + d = 0 (a不为0)" << std::endl;
				std::cout << "请输入系数a, b, c, d：(无需输入Ctrl+Z以确认)" << std::endl;
				std::cin >> a >> b >> c >> d;
				Linear::Equation::Cubic(a, b, c, d, result);
				Linear::Equation::Print(result);
			}
			static void CT3()
			{

			}
			static void CT4()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::vector<double> result;
					std::cout << "请输入未知数数量：" << std::endl;
					std::cin >> n;
					std::cout << "请输入系数方阵：" << std::endl;
					coeff.Input(n);
					std::cout << "请输入右侧值：" << std::endl;
					Data::Process::Input(right_value, n);
					Linear::Equation::Set(coeff, right_value, result);
					Linear::Equation::Print(result);
				}
				if (YN == 'Y' || YN == 'y') {
					std::vector<double> result;
					std::cout << "已继承上次输入的系数矩阵。" << std::endl;
					std::cout << "请输入右侧值：" << std::endl;
					Data::Process::Input(right_value, coeff.Rank());
					Linear::Equation::Set(coeff, right_value, result);
					Linear::Equation::Print(result);
				}
			}

		};

		class CS2
		{
		public:
			static void CT1()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入向量维数：" << std::endl;
					std::cin >> n;
					std::cout << "请输入向量x：" << std::endl;
					x.Input(n);
					std::cout << "向量x的长度为：" << std::endl << x.Length() << std::endl;
				}
				if (YN == 'Y' || YN == 'y') {
					std::cout << "向量x的长度为：" << std::endl << x.Length() << std::endl;
				}
			}
			static void CT2()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入向量维数：" << std::endl;
					std::cin >> n;
					std::cout << "请输入向量x：" << std::endl;
					x.Input(n);
					std::cout << "请输入向量y：" << std::endl;
					y.Input(n);
					std::cout << "向量x, y的标量积为：" << std::endl << Linear::Vector::ScalarPro(x, y) << std::endl;
				}
				if (YN == 'Y' || YN == 'y') {
					std::cout << "向量x, y的标量积为：" << std::endl << Linear::Vector::ScalarPro(x, y) << std::endl;
				}
			}
			static void CT3()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "向量维数固定为3，请输入向量x：" << std::endl;
					x.Input(3);
					std::cout << "请输入向量y：" << std::endl;
					y.Input(3);
					Linear::Vector result = x.VectorPro(y);
					std::cout << "向量x, y的向量积为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Vector result = x.VectorPro(y);
					std::cout << "向量x, y的向量积为：" << std::endl;
					result.Print();
				}
			}
			static void CT4()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入向量维数：" << std::endl;
					std::cin >> n;
					std::cout << "请输入向量x：" << std::endl;
					x.Input(n);
					std::cout << "请输入向量y：" << std::endl;
					y.Input(n);
					std::cout << "向量x, y的夹角为：" << std::endl << Linear::Vector::Angle(x, y) << std::endl;
				}
				if (YN == 'Y' || YN == 'y') {
					std::cout << "向量x, y的夹角为：" << std::endl << Linear::Vector::Angle(x, y) << std::endl;
				}
			}

		};

		class CS3
		{
		public:
			static void CT1()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入行列式的阶：" << std::endl;
					std::cin >> n;
					std::cout << "请输入行列式：" << std::endl;
					matrix.Input(n);
					std::cout << "行列式的值为：" << std::endl << Linear::Matrix::Det(matrix) << std::endl;
				}
				if (YN == 'Y' || YN == 'y') {
					std::cout << "行列式的值为：" << std::endl << Linear::Matrix::Det(matrix) << std::endl;
				}
			}
			static void CT2()
			{
				if (YN == 'N' || YN == 'n') {
					int R = 0, C = 0;
					std::cout << "请输入矩阵的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵1：" << std::endl;
					matrix1.Input(R, C);
					std::cout << "请输入矩阵2：" << std::endl;
					matrix2.Input(R, C);
					Linear::Matrix result = matrix1 + matrix2;
					std::cout << "和矩阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix1 + matrix2;
					std::cout << "和矩阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT3()
			{
				if (YN == 'N' || YN == 'n') {
					int R = 0, C = 0;
					std::cout << "请输入矩阵的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵1：" << std::endl;
					matrix1.Input(R, C);
					std::cout << "请输入矩阵2：" << std::endl;
					matrix2.Input(R, C);
					Linear::Matrix result = matrix1 - matrix2;
					std::cout << "差矩阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix1 - matrix2;
					std::cout << "差矩阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT4()
			{
				if (YN == 'N' || YN == 'n') {
					int R = 0, C = 0;
					std::cout << "请输入矩阵1的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵1的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵1：" << std::endl;
					matrix1.Input(R, C);
					std::cout << "请输入矩阵2的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵2的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵2：" << std::endl;
					matrix2.Input(R, C);
					Linear::Matrix result = matrix1 * matrix2;
					std::cout << "矩阵之积为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix1 * matrix2;
					std::cout << "矩阵之积为：" << std::endl;
					result.Print();
				}
			}
			static void CT5()
			{
				if (YN == 'N' || YN == 'n') {
					int R = 0, C = 0;
					std::cout << "请输入矩阵的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵：" << std::endl;
					matrix1.Input(R, C);
					Linear::Matrix result = matrix1.GausElmn();
					std::cout << "矩阵的行阶梯形矩阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix1.GausElmn();
					std::cout << "矩阵的行阶梯形矩阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT6()
			{
				if (YN == 'N' || YN == 'n') {
					int R = 0, C = 0;
					std::cout << "请输入矩阵的行数R：" << std::endl;
					std::cin >> R;
					std::cout << "请输入矩阵的列数C：" << std::endl;
					std::cin >> C;
					std::cout << "请输入矩阵：" << std::endl;
					matrix1.Input(R, C);
					Linear::Matrix result = matrix1.RowSim();
					std::cout << "矩阵的行最简形矩阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix1.RowSim();
					std::cout << "矩阵的行最简形矩阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT7()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0, power = 0;
					std::cout << "请输入方阵的阶：" << std::endl;
					std::cin >> n;
					std::cout << "请输入方阵：" << std::endl;
					matrix.Input(n);
					std::cout << "请输入方阵的幂：" << std::endl;
					std::cin >> power;
					Linear::Matrix result = matrix.Power(power);
					std::cout << "方阵的" << power << "次幂为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					int power = 0;
					std::cout << "已继承上次输入的方阵。" << std::endl;
					std::cout << "请输入方阵的幂：" << std::endl;
					std::cin >> power;
					Linear::Matrix result = matrix.Power(power);
					std::cout << "方阵的" << power << "次幂为：" << std::endl;
					result.Print();
				}
			}
			static void CT8()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入方阵的阶：" << std::endl;
					std::cin >> n;
					std::cout << "请输入方阵：" << std::endl;
					matrix.Input(n);
					Linear::Matrix result = matrix.Adjoint();
					std::cout << "方阵的伴随阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix.Adjoint();
					std::cout << "方阵的伴随阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT9()
			{
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入方阵的阶：" << std::endl;
					std::cin >> n;
					std::cout << "请输入方阵：" << std::endl;
					matrix.Input(n);
					Linear::Matrix result = matrix.Inverse();
					std::cout << "方阵的逆阵为：" << std::endl;
					result.Print();
				}
				if (YN == 'Y' || YN == 'y') {
					Linear::Matrix result = matrix.Inverse();
					std::cout << "方阵的逆阵为：" << std::endl;
					result.Print();
				}
			}
			static void CT10()
			{
				std::vector<double> result;
				result.push_back(0);
				if (YN == 'N' || YN == 'n') {
					int n = 0;
					std::cout << "请输入方阵的阶：" << std::endl;
					std::cin >> n;
					std::cout << "请输入方阵：" << std::endl;
					matrix.Input(n);
					for (int i = -100; i <= 100; i++) {
						Linear::Matrix temp = matrix;
						Linear::Matrix unit = unit.UnitMatrix(n, i);
						temp = temp - unit;
						if (fabs(Linear::Matrix::Det(temp)) <= 1e-6) result.push_back(i);
					}
					std::cout << "方阵的特征值为：" << std::endl;
					if (result.size() == 1) std::cout << "无+-100以内整数值。" << std::endl;
					else Linear::Equation::Print(result);
				}
				if (YN == 'Y' || YN == 'y') {
					for (int i = -100; i <= 100; i++) {
						Linear::Matrix temp = matrix;
						Linear::Matrix unit = unit.UnitMatrix(matrix.Rank(), i);
						temp = temp - unit;
						if (fabs(Linear::Matrix::Det(temp)) <= 1e-6) result.push_back(i);
					}
					std::cout << "方阵的特征值为：" << std::endl;
					if (result.size() == 1) std::cout << "无+-100以内整数值。" << std::endl;
					else Linear::Equation::Print(result);
				}
			}
			static void CT11()
			{

			}

		};

	};

};