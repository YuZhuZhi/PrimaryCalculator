#pragma once
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>
#define Max(a, b) (a > b) ? a : b;
#define Min(a, b) (a > b) ? b : a;

class Data
{
public:
	class Process //数据处理
	{
	public:
		static void Input(std::vector<double>& data_set) //数据集输入：以vector为载体，为了契合Matrix，0位置不填充有效数据
		{
			data_set.clear();
			data_set.push_back(0); //0位置以0填充
			double temp = 0;
			while (std::cin >> temp) { //当未输入Ctrl+Z前，总是接受输入
				data_set.push_back(temp);
			}
			std::cin.clear(); //清除Ctrl+Z状态
		}
		static void Input(std::vector<double>& data_set, int times) //数据集输入重载：给定数据容量，不需以Ctrl+Z结束输入
		{
			data_set.clear();
			data_set.push_back(0); //0位置以0填充
			double temp = 0;
			for (int i = 1; i <= times; i++) {
				std::cin >> temp;
				data_set.push_back(temp);
			}
		}
		static void Print(std::vector<double> data_set) //数据集输出：中间以空格间隔
		{
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				std::cout << *vdit << " ";
			}
			std::cout << std::endl << std::endl;
		}
		static int Capacity(std::vector<double> data_set) //返回数据容量
		{
			return (data_set.size() - 1); //0位置是无效数据
		}
		static double Summary(std::vector<double> data_set) //计算数据集总和
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + *vdit;
			}
			return sum;
		}
		static double Summary(std::vector<double> data_set, int power) //计算数据集中元素power次幂后的总和
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + pow(*vdit, power);
			}
			return sum;
		}
		static double CroseSum(std::vector<double> x, std::vector<double> y) //计算两个等容数据集对应位置元素之积的总和，称为交叉和
		{
			if (x.size() != y.size()) return (0);
			double crssum = 0;
			for (int i = 1; i <= x.size() - 1; i++) {
				crssum = crssum + x[i] * y[i];
			}
			return (crssum);
		}
		static double Average(std::vector<double> data_set) //计算数据集的平均数
		{
			double sum = Summary(data_set);
			return (sum / (double)(Capacity(data_set)));
		}
		static double CroseAver(std::vector<double> data_set1, std::vector<double> data_set2) //计算两个等容数据集的交叉平均数
		{
			return (CroseSum(data_set1, data_set2) / (double)(Capacity(data_set)));
		}
		static double SquareAver(std::vector<double> data_set) //计算数据集中元素平方后的平均数
		{
			return (Summary(data_set, 2) / (double)(Capacity(data_set)));
		}
		static double StdDev(std::vector<double> data_set) //计算数据集的标准差
		{
			double aver = Average(data_set);
			double dif_sq = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				dif_sq = dif_sq + (aver - *vdit) * (aver - *vdit);
			}
			return (sqrt(dif_sq / (double)(Capacity(data_set))));
		}
		static double AverStdDev(std::vector<double> data_set) //计算数据集的平均值标准差
		{
			return (StdDev(data_set) / sqrt((double)(Capacity(data_set) - 1)));
		}
		static double MedianNum(std::vector<double> data_set) //计算数据集的中位数
		{
			sort(data_set.begin() + 1, data_set.end()); //排序，并且忽略0位置
			int size = Capacity(data_set); //计算数据容量后分奇偶
			if (size % 2 == 0) {
				return ((data_set[(size / 2)] + data_set[(size / 2) + 1]) / 2.0);
			}
			else return (data_set[(size / 2) + 1]);
		}
		static int ModeNum(std::vector<double> data_set, std::vector<int>& mode_nums) //计算数据集的众数，用mode_nums承装找到的众数，返回值是其中一个众数出现的次数
		{
			mode_nums.resize(0);
			data_set.erase(data_set.begin()); //删除0位置的0元素
			std::vector<int> data_set_int(data_set.begin(), data_set.end()); //如果输入的数据集是浮点数，转化为整型
			sort(data_set_int.begin(), data_set_int.end()); //排序
			int max = 0, count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //如果当前数与前一个数相等，那么计数器count+1
				else {
					max = Max(max, count); //否则比较max和count并更新max为较大者
					count = 1;
				}
				if (viit == data_set_int.end() - 1) max = Max(max, count); //当当前数是最后一个数，还需要再比较一次
			}
			count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //如果当前数与前一个数相等，那么计数器count+1
				else {
					if (count == max) mode_nums.push_back(*(viit - 1)); //如果当前数的个数恰是max，那么当前数就是众数之一
					count = 1;
				}
				if (viit == data_set_int.end() - 1 && count == max) mode_nums.push_back(*(viit - 1)); //当当前数是最后一个数，还需要再比较一次
			}
			return max;
		}
	};

	class LinearFit //线性拟合
	{
	public:
		static void Input(std::vector<double>& x) //自/因变量输入
		{
			Process::Input(x); //衍生自数据集输入
		}
		static void Input(std::vector<double>& x, int times) //指定数据容量的自/因变量输入
		{
			Process::Input(x, times);
		}
		static void Print(std::vector<double> x) //自/因变量输出
		{
			Process::Print(x);
		}
		static double K(std::vector<double> x, std::vector<double> y) //计算拟合直线的斜率
		{
			int n = x.size() - 1;
			if (y.size() - 1 != n) return 0;
			else return ((Process::CroseSum(x, y) - (double)n * Process::Average(x) * Process::Average(y)) / (Process::Summary(x, 2) - (double)n * Process::Average(x) * Process::Average(x)));
		}
		static double B(std::vector<double> x, std::vector<double> y) //计算拟合直线的截距
		{
			return (Process::Average(y) - K(x, y) * Process::Average(x));
		}
		static double R(std::vector<double> x, std::vector<double> y) //计算拟合度
		{
			double x_sqaver = Process::SquareAver(x), y_sqaver = Process::SquareAver(y);
			double x_aver = Process::Average(x), y_aver = Process::Average(y);
			return ((Process::CroseAver(x, y) - x_aver * y_aver) / (sqrt((x_sqaver - x_aver * x_aver) * (y_sqaver - y_aver * y_aver))));
		}
	};

};