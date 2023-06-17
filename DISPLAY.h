#pragma once
#include <iostream>
#include <cstdio>
#include <Windows.h>
#include "OPERATE.h"

const int sleeptime = 200;

extern char YN;
extern int choicefir, choicesec, choicethr;
extern std::vector<double> data_set, x, y, right_value;
extern std::vector<std::vector<double> > matrix, matrix1, matrix2, matrix3, coeff;

class Display
{
public:
	static void Guide()
	{
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "欢迎使用数学计算器ver0.56(TEST)。(Programmed by 御伫之)" << std::endl;
		std::cout << "(御伫之 知乎主页：https://www.zhihu.com/people/yuzhuzhi)" << std::endl;
		std::cout << "(本计算器开源于：https://github.com/YuZhuZhi/PrimaryCalculator)" << std::endl;
		std::cout << "使用方法：根据提示输入数字或者字母后回车以选择操作(字母可不区分大小写)。" << std::endl;
		std::cout << "输入数据时可以用空格间隔或者每输入一次即按一次回车以提交。" << std::endl;
		std::cout << "也可以从Excel中复制数据后提交。" << std::endl;
		std::cout << "结束输入方法：提交完最后一个数据后，按Ctrl+Z出现\"^Z\"后提交即可。示例：" << std::endl;
		std::cout << "1 2 3 4 5" << std::endl << "^Z" << std::endl;
		std::cout << "继承上次输入是指继续使用上次已输入数据，再带到另一种或者同一种操作中。" << std::endl;
		std::cout << "1e5表示1乘10的5次方，以此类推。" << std::endl;
		std::cout << "祝您使用愉快。" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
	}
	static void ChooseYN()
	{
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "是否继承上一次输入的数据？" << std::endl;
		std::cout << "Y\t是" << std::endl;
		std::cout << "N\t否" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> YN;
	}
	static void ChooseFir()
	{
		system("color F0");
		//Sleep(sleeptime);
		if (YN == 'N' || YN == 'n') ChooseFirWithN();
		if (YN == 'Y' || YN == 'y') ChooseFirWithY();
		ChooseYN();
	}
	static void ChooseFirWithN()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：" << std::endl;
		std::cout << "1\t常规数据处理(6)" << std::endl;
		std::cout << "2\t数学运算与多项式(5)(暂未搭载)" << std::endl;
		std::cout << "3\t线性代数(3)" << std::endl << std::endl;
		std::cout << "4\t清除已显示内容" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicefir;
		if (std::cin.fail()) {
			std::cin.clear();
			std::cin.ignore();
			std::cout << "搞神马呢小伙计?" << std::endl;
			std::cout << "请重新选择：" << std::endl;
			ChooseFir();
		}
		switch (choicefir)
		{
		case 1:
			ChooseSec1();
			break;
		case 2:
			ChooseSec2();
			break;
		case 3:
			ChooseSec3();
			break;
		case 4:
			system("cls");
			ChooseFirWithN();
			break;
		default:
			ChooseFir();
			break;
		}
	}
	static void ChooseFirWithY()
	{
		switch (choicefir)
		{
		case 1:
			ChooseSec1();
			break;
		case 2:
			ChooseSec2();
			break;
		case 3:
			ChooseSec3();
			break;
		default:
			ChooseFirWithY();
			break;
		}
	}
	static void ChooseSec1()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：" << std::endl;
		std::cout << "常规数据处理：" << std::endl;
		std::cout << "1\t求和" << std::endl;
		std::cout << "2\t求均" << std::endl;
		std::cout << "3\t求标准差(数) / 实验标准差(物)" << std::endl;
		std::cout << "4\t求平均值标准差(数) / 平均值的实验标准差(物)" << std::endl;
		std::cout << "5\t求直接测量量的不确定度" << std::endl;
		std::cout << "6\t单变量线性拟合" << std::endl << std::endl;
		std::cout << "7...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicesec;
		switch (choicesec)
		{
		case 1:
			Operate::CF1::CS1();
			break;
		case 2:
			Operate::CF1::CS2();
			break;
		case 3:
			Operate::CF1::CS3();
			break;
		case 4:
			Operate::CF1::CS4();
			break;
		case 5:
			Operate::CF1::CS5();
			break;
		case 6:
			Operate::CF1::CS6();
			break;
		default:
			ChooseFir();
			break;
		}
	}
	static void ChooseSec2()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：" << std::endl;
		std::cout << "数学运算与多项式：" << std::endl;
		std::cout << "1\t阶乘" << std::endl;
		std::cout << "2\t双阶乘" << std::endl;
		std::cout << "3\t36以内进制转换" << std::endl;
		std::cout << "4\t已展开二多项式之积" << std::endl;
		std::cout << "5\t可完全分解的多项式的展开" << std::endl << std::endl;
		std::cout << "6...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicesec;
		switch (choicesec)
		{
		case 1:
			Operate::CF2::CS1();
			break;
		case 2:
			Operate::CF2::CS2();
			break;
		case 3:
			Operate::CF2::CS3();
			break;
		case 4:
			Operate::CF2::CS4();
			break;
		case 5:
			Operate::CF2::CS5();
			break;
		default:
			ChooseFir();
			break;
		}
	}
	static void ChooseSec3()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：(输入数据时无需输入Ctrl+Z以确认)" << std::endl;
		std::cout << "线性代数:" << std::endl;
		std::cout << "1\t解方程(组)(4)" << std::endl;
		std::cout << "2\t向量(4)" << std::endl;
		std::cout << "3\t行列式与矩阵(10)" << std::endl << std::endl;
		std::cout << "4...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicesec;
		switch (choicesec)
		{
		case 1:
			ChooseThr1();
			break;
		case 2:
			ChooseThr2();
			break;
		case 3:
			ChooseThr3();
			break;
		default:
			ChooseFir();
			break;
		}
	}
	static void ChooseThr1()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：(输入时无需输入Ctrl+Z以确认)" << std::endl;
		std::cout << "解方程(组):" << std::endl;
		std::cout << "1\t解一元二次方程" << std::endl;
		std::cout << "2\t解一元三次方程(暂未搭载)" << std::endl;
		std::cout << "3\t解一元n次方程(暂未搭载)" << std::endl;
		std::cout << "4\t解n元方程组" << std::endl << std::endl;
		std::cout << "5...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicethr;
		switch (choicethr)
		{
		case 1:
			Operate::CF3::CS1::CT1();
			break;
		case 2:
			Operate::CF3::CS1::CT2();
			break;
		case 3:
			Operate::CF3::CS1::CT3();
			break;
		case 4:
			Operate::CF3::CS1::CT4();
			break;
		default:
			ChooseSec3();
			break;
		}
	}
	static void ChooseThr2()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：(输入时无需输入Ctrl+Z以确认)" << std::endl;
		std::cout << "向量：" << std::endl;
		std::cout << "1\t求向量的长度" << std::endl;
		std::cout << "2\t求二向量的标量积" << std::endl;
		std::cout << "3\t求二向量的向量积" << std::endl;
		std::cout << "4\t求二向量的夹角(返回弧度)" << std::endl << std::endl;
		std::cout << "5...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicethr;
		switch (choicethr)
		{
		case 1:
			Operate::CF3::CS2::CT1();
			break;
		case 2:
			Operate::CF3::CS2::CT2();
			break;
		case 3:
			Operate::CF3::CS2::CT3();
			break;
		case 4:
			Operate::CF3::CS2::CT4();
			break;
		default:
			ChooseSec3();
			break;
		}
	}
	static void ChooseThr3()
	{
		Sleep(sleeptime);
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "选择你需要的操作：(输入时无需输入Ctrl+Z以确认)" << std::endl;
		std::cout << "(行列式/方阵可互相继承，但不可继承到矩阵；反之亦然)" << std::endl;
		std::cout << "行列式与矩阵:" << std::endl;
		std::cout << "1\t求行列式的值" << std::endl;
		std::cout << "2\t求二矩阵之和" << std::endl;
		std::cout << "3\t求二矩阵之差" << std::endl;
		std::cout << "4\t求二矩阵之积" << std::endl;
		std::cout << "5\t求矩阵的行阶梯形矩阵(高斯消元法)" << std::endl;
		std::cout << "6\t求矩阵的行最简形矩阵" << std::endl;
		std::cout << "7\t求方阵的n次幂" << std::endl;
		std::cout << "8\t求方阵的伴随阵" << std::endl;
		std::cout << "9\t求方阵的逆阵" << std::endl;
		std::cout << "10\t求方阵的特征值(初步)" << std::endl;
		std::cout << "11\t方阵施密特正交化(暂未搭载)" << std::endl << std::endl;
		std::cout << "12...\t回到上一步" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicethr;
		switch (choicethr)
		{
		case 1:
			Operate::CF3::CS3::CT1();
			break;
		case 2:
			Operate::CF3::CS3::CT2();
			break;
		case 3:
			Operate::CF3::CS3::CT3();
			break;
		case 4:
			Operate::CF3::CS3::CT4();
			break;
		case 5:
			Operate::CF3::CS3::CT5();
			break;
		case 6:
			Operate::CF3::CS3::CT6();
			break;
		case 7:
			Operate::CF3::CS3::CT7();
			break;
		case 8:
			Operate::CF3::CS3::CT8();
			break;
		case 9:
			Operate::CF3::CS3::CT9();
			break;
		case 10:
			Operate::CF3::CS3::CT10();
			break;
		case 11:
			Operate::CF3::CS3::CT11();
			break;
		default:
			ChooseSec3();
			break;
		}
	}

};