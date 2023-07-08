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
	class Process //���ݴ���
	{
	public:
		static void Input(std::vector<double>& data_set) //���ݼ����룺��vectorΪ���壬Ϊ������Matrix��0λ�ò������Ч����
		{
			data_set.clear();
			data_set.push_back(0); //0λ����0���
			double temp = 0;
			while (std::cin >> temp) { //��δ����Ctrl+Zǰ�����ǽ�������
				data_set.push_back(temp);
			}
			std::cin.clear(); //���Ctrl+Z״̬
		}
		static void Input(std::vector<double>& data_set, int times) //���ݼ��������أ���������������������Ctrl+Z��������
		{
			data_set.clear();
			data_set.push_back(0); //0λ����0���
			double temp = 0;
			for (int i = 1; i <= times; i++) {
				std::cin >> temp;
				data_set.push_back(temp);
			}
		}
		static void Print(std::vector<double> data_set) //���ݼ�������м��Կո���
		{
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				std::cout << *vdit << " ";
			}
			std::cout << std::endl << std::endl;
		}
		static int Capacity(std::vector<double> data_set) //������������
		{
			return (data_set.size() - 1); //0λ������Ч����
		}
		static double Summary(std::vector<double> data_set) //�������ݼ��ܺ�
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + *vdit;
			}
			return sum;
		}
		static double Summary(std::vector<double> data_set, int power) //�������ݼ���Ԫ��power���ݺ���ܺ�
		{
			double sum = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				sum = sum + pow(*vdit, power);
			}
			return sum;
		}
		static double CroseSum(std::vector<double> x, std::vector<double> y) //���������������ݼ���Ӧλ��Ԫ��֮�����ܺͣ���Ϊ�����
		{
			if (x.size() != y.size()) return (0);
			double crssum = 0;
			for (int i = 1; i <= x.size() - 1; i++) {
				crssum = crssum + x[i] * y[i];
			}
			return (crssum);
		}
		static double Average(std::vector<double> data_set) //�������ݼ���ƽ����
		{
			double sum = Summary(data_set);
			return (sum / (double)(Capacity(data_set)));
		}
		static double CroseAver(std::vector<double> data_set1, std::vector<double> data_set2) //���������������ݼ��Ľ���ƽ����
		{
			return (CroseSum(data_set1, data_set2) / (double)(Capacity(data_set1)));
		}
		static double SquareAver(std::vector<double> data_set) //�������ݼ���Ԫ��ƽ�����ƽ����
		{
			return (Summary(data_set, 2) / (double)(Capacity(data_set)));
		}
		static double StdDev(std::vector<double> data_set) //�������ݼ��ı�׼��
		{
			double aver = Average(data_set);
			double dif_sq = 0;
			for (std::vector<double>::iterator vdit = data_set.begin() + 1; vdit != data_set.end(); vdit++) {
				dif_sq = dif_sq + (aver - *vdit) * (aver - *vdit);
			}
			return (sqrt(dif_sq / (double)(Capacity(data_set))));
		}
		static double AverStdDev(std::vector<double> data_set) //�������ݼ���ƽ��ֵ��׼��
		{
			return (StdDev(data_set) / sqrt((double)(Capacity(data_set) - 1)));
		}
		static double MedianNum(std::vector<double> data_set) //�������ݼ�����λ��
		{
			sort(data_set.begin() + 1, data_set.end()); //���򣬲��Һ���0λ��
			int size = Capacity(data_set); //�����������������ż
			if (size % 2 == 0) {
				return ((data_set[(size / 2)] + data_set[(size / 2) + 1]) / 2.0);
			}
			else return (data_set[(size / 2) + 1]);
		}
		static int ModeNum(std::vector<double> data_set, std::vector<int>& mode_nums) //�������ݼ�����������mode_nums��װ�ҵ�������������ֵ������һ���������ֵĴ���
		{
			mode_nums.resize(0);
			data_set.erase(data_set.begin()); //ɾ��0λ�õ�0Ԫ��
			std::vector<int> data_set_int(data_set.begin(), data_set.end()); //�����������ݼ��Ǹ�������ת��Ϊ����
			sort(data_set_int.begin(), data_set_int.end()); //����
			int max = 0, count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //�����ǰ����ǰһ������ȣ���ô������count+1
				else {
					max = Max(max, count); //����Ƚ�max��count������maxΪ�ϴ���
					count = 1;
				}
				if (viit == data_set_int.end() - 1) max = Max(max, count); //����ǰ�������һ����������Ҫ�ٱȽ�һ��
			}
			count = 1;
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //�����ǰ����ǰһ������ȣ���ô������count+1
				else {
					if (count == max) mode_nums.push_back(*(viit - 1)); //�����ǰ���ĸ���ǡ��max����ô��ǰ����������֮һ
					count = 1;
				}
				if (viit == data_set_int.end() - 1 && count == max) mode_nums.push_back(*(viit - 1)); //����ǰ�������һ����������Ҫ�ٱȽ�һ��
			}
			return max;
		}
	};

	class LinearFit //�������
	{
	public:
		static void Input(std::vector<double>& x) //��/���������
		{
			Process::Input(x); //���������ݼ�����
		}
		static void Input(std::vector<double>& x, int times) //ָ��������������/���������
		{
			Process::Input(x, times);
		}
		static void Print(std::vector<double> x) //��/��������
		{
			Process::Print(x);
		}
		static double K(std::vector<double> x, std::vector<double> y) //�������ֱ�ߵ�б��
		{
			int n = x.size() - 1;
			if (y.size() - 1 != n) return 0;
			else return ((Process::CroseSum(x, y) - (double)n * Process::Average(x) * Process::Average(y)) / (Process::Summary(x, 2) - (double)n * Process::Average(x) * Process::Average(x)));
		}
		static double B(std::vector<double> x, std::vector<double> y) //�������ֱ�ߵĽؾ�
		{
			return (Process::Average(y) - K(x, y) * Process::Average(x));
		}
		static double R(std::vector<double> x, std::vector<double> y) //������϶�
		{
			double x_sqaver = Process::SquareAver(x), y_sqaver = Process::SquareAver(y);
			double x_aver = Process::Average(x), y_aver = Process::Average(y);
			return ((Process::CroseAver(x, y) - x_aver * y_aver) / (sqrt((x_sqaver - x_aver * x_aver) * (y_sqaver - y_aver * y_aver))));
		}
	};

};