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
	class Process
	{
	public:
		static void Input(std::vector<double>& data_set)
		{
			data_set.clear();
			data_set.push_back(0);
			double temp = 0;
			while (std::cin >> temp) {
				data_set.push_back(temp);
			}
			std::cin.clear();
		}
		static void Input(std::vector<double>& data_set, int times)
		{
			data_set.clear();
			data_set.push_back(0);
			double temp = 0;
			for (int i = 1; i <= times; i++) {
				std::cin >> temp;
				data_set.push_back(temp);
			}
		}
		static void Print(std::vector<double> data_set)
		{
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				std::cout << *vdit << " ";
			}
			std::cout << std::endl << std::endl;
		}
		static int Capacity(std::vector<double> data_set)
		{
			return (data_set.size() - 1);
		}
		static double Summary(std::vector<double> data_set)
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + *vdit;
			}
			return sum;
		}
		static double Summary(std::vector<double> data_set, int power)
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + pow(*vdit, power);
			}
			return sum;
		}
		static double CroseSum(std::vector<double> x, std::vector<double> y)
		{
			if (x.size() != y.size()) return (0);
			double crssum = 0;
			for (int i = 1; i <= x.size() - 1; i++) {
				crssum = crssum + x[i] * y[i];
			}
			return (crssum);
		}
		static double Average(std::vector<double> data_set)
		{
			double sum = Summary(data_set);
			return (sum / (double)(data_set.size() - 1));
		}
		static double CroseAver(std::vector<double> data_set1, std::vector<double> data_set2)
		{
			return (CroseSum(data_set1, data_set2) / (double)(data_set1.size() - 1));
		}
		static double SquareAver(std::vector<double> data_set)
		{
			return (Summary(data_set, 2) / (double)(data_set.size() - 1));
		}
		static double StdDev(std::vector<double> data_set)
		{
			double aver = Average(data_set);
			double dif_sq = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				dif_sq = dif_sq + (aver - *vdit) * (aver - *vdit);
			}
			return (sqrt(dif_sq / (double)(data_set.size() - 1)));
		}
		static double AverStdDev(std::vector<double> data_set)
		{
			return (StdDev(data_set) / sqrt((double)(data_set.size() - 2)));
		}
		static double MedianNum(std::vector<double> data_set)
		{
			sort(data_set.begin() + 1, data_set.end());
			int size = data_set.size() - 1;
			if (size % 2 == 0) {
				return ((data_set[(size / 2)] + data_set[(size / 2) + 1]) / 2.0);
			}
			else return (data_set[(size / 2) + 1]);
		}
		static int ModeNum(std::vector<double> data_set, std::vector<int>& mode_nums)
		{
			mode_nums.resize(0);
			data_set.erase(data_set.begin());
			std::vector<int> data_set_int(data_set.begin(), data_set.end());
			sort(data_set_int.begin(), data_set_int.end());
			int max = 0, count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++;
				else {
					max = Max(max, count);
					count = 1;
				}
				if (viit == data_set_int.end() - 1) max = Max(max, count);
			}
			count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++;
				else {
					if (count == max) mode_nums.push_back(*(viit - 1));
					count = 1;
				}
				if (viit == data_set_int.end() - 1 && count == max) mode_nums.push_back(*(viit - 1));
			}
			return max;
		}
	};

	class LinearFit
	{
	public:
		static void Input(std::vector<double>& x)
		{
			Process::Input(x);
		}
		static void Input(std::vector<double>& x, int times)
		{
			Process::Input(x, times);
		}
		static void Print(std::vector<double> x)
		{
			Process::Print(x);
		}
		static double K(std::vector<double> x, std::vector<double> y)
		{
			int n = x.size() - 1;
			if (y.size() - 1 != n) return 0;
			else return ((Process::CroseSum(x, y) - (double)n * Process::Average(x) * Process::Average(y)) / (Process::Summary(x, 2) - (double)n * Process::Average(x) * Process::Average(x)));
		}
		static double B(std::vector<double> x, std::vector<double> y)
		{
			return (Process::Average(y) - K(x, y) * Process::Average(x));
		}
		static double R(std::vector<double> x, std::vector<double> y)
		{
			double x_sqaver = Process::SquareAver(x), y_sqaver = Process::SquareAver(y);
			double x_aver = Process::Average(x), y_aver = Process::Average(y);
			return ((Process::CroseAver(x, y) - x_aver * y_aver) / (sqrt((x_sqaver - x_aver * x_aver) * (y_sqaver - y_aver * y_aver))));
		}
	};

};