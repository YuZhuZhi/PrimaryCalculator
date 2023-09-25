#pragma once
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

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
		static void Print(const std::vector<double>& data_set) //数据集输出：中间以空格间隔
		{
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				std::cout << *vdit << " ";
			}
			std::cout << std::endl << std::endl;
		}
		static int Capacity(const std::vector<double>& data_set) //返回数据容量
		{
			return (data_set.size() - 1); //0位置是无效数据
		}
		static double Summary(const std::vector<double>& data_set) //计算数据集总和
		{
			double sum = 0;
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				sum = sum + *vdit;
			}
			return sum;
		}
		static double Summary(const std::vector<double>& data_set, const int power) //计算数据集中元素power次幂后的总和
		{
			double sum = 0;
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				sum = sum + pow(*vdit, power);
			}
			return sum;
		}
		static double CroseSum(const std::vector<double>& x, const std::vector<double>& y) //计算两个等容数据集对应位置元素之积的总和，称为交叉和
		{
			if (Capacity(x) != Capacity(y)) return (0);
			double crssum = 0;
			for (int i = 1; i <= Capacity(x); i++) {
				crssum = crssum + x[i] * y[i];
			}
			return (crssum);
		}
		static double Average(const std::vector<double>& data_set) //计算数据集的平均数
		{
			double sum = Summary(data_set);
			return (sum / (double)(Capacity(data_set)));
		}
		static double CroseAver(const std::vector<double>& data_set1, const std::vector<double>& data_set2) //计算两个等容数据集的交叉平均数
		{
			return (CroseSum(data_set1, data_set2) / (double)(Capacity(data_set1)));
		}
		static double SquareAver(const std::vector<double>& data_set) //计算数据集中元素平方后的平均数
		{
			return (Summary(data_set, 2) / (double)(Capacity(data_set)));
		}
		static double StdDev(const std::vector<double>& data_set) //计算数据集的标准差
		{
			std::vector<double> temp = data_set;
			double aver = Average(temp);
			double dif_sq = 0;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				dif_sq = dif_sq + (aver - *vdit) * (aver - *vdit);
			}
			return (sqrt(dif_sq / (double)(Capacity(data_set))));
		}
		static double AverStdDev(const std::vector<double>& data_set) //计算数据集的平均值标准差
		{
			return (StdDev(data_set) / sqrt((double)(Capacity(data_set) - 1)));
		}
		static double MedianNum(const std::vector<double>& data_set) //计算数据集的中位数
		{
			std::vector<double> temp = data_set;
			sort(temp.begin() + 1, temp.end()); //排序，并且忽略0位置
			int size = Capacity(temp); //计算数据容量后分奇偶
			if (size % 2 == 0) {
				return ((temp[(size / 2)] + temp[(size / 2) + 1]) / 2.0);
			}
			else return (temp[(size / 2) + 1]);
		}
		static int ModeNum(const std::vector<double>& data_set, std::vector<int>& mode_nums) //计算数据集的众数，用mode_nums承装找到的众数，返回值是其中一个众数出现的次数
		{
			int count = 1;
			struct dat {
				int _num;
				int _count;
			};
			std::vector<dat> stat;
			mode_nums.resize(0);
			std::vector<int> data_set_int(data_set.begin() + 1, data_set.end()); //如果输入的数据集是浮点数，转化为整型
			sort(data_set_int.begin(), data_set_int.end()); //排序
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //如果当前数与前一个数相等，那么计数器count+1
				else {
					stat.push_back({ *(viit - 1), count });
					count = 1; //计数器归1
				}
				if (viit == data_set_int.end() - 1) stat.push_back({ *(viit - 1), count }); //当当前数是最后一个数，还需要再push一次
			}
			sort(stat.begin(), stat.end(), [](dat x, dat y) {return x._count > y._count; }); //按dat中的_count大小降序排列
			const int max = stat[0]._count;
			for (int i = 0; stat[i]._count == max; i++) mode_nums.push_back(stat[i]._num);
			return stat[0]._count;
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
		static void Print(const std::vector<double>& x) //自/因变量输出
		{
			Process::Print(x);
		}
		static double K(const std::vector<double>& x, const std::vector<double>& y) //计算拟合直线的斜率
		{
			int n = x.size() - 1;
			if (y.size() - 1 != n) return 0;
			else return ((Process::CroseSum(x, y) - (double)n * Process::Average(x) * Process::Average(y)) / (Process::Summary(x, 2) - (double)n * Process::Average(x) * Process::Average(x)));
		}
		static double B(const std::vector<double>& x, const std::vector<double>& y) //计算拟合直线的截距
		{
			return (Process::Average(y) - K(x, y) * Process::Average(x));
		}
		static double R(const std::vector<double>& x, const std::vector<double>& y) //计算拟合度
		{
			double x_sqaver = Process::SquareAver(x), y_sqaver = Process::SquareAver(y);
			double x_aver = Process::Average(x), y_aver = Process::Average(y);
			return ((Process::CroseAver(x, y) - x_aver * y_aver) / (sqrt((x_sqaver - x_aver * x_aver) * (y_sqaver - y_aver * y_aver))));
		}
	};

};