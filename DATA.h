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
		static void Print(const std::vector<double>& data_set) //���ݼ�������м��Կո���
		{
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				std::cout << *vdit << " ";
			}
			std::cout << std::endl << std::endl;
		}
		static int Capacity(const std::vector<double>& data_set) //������������
		{
			return (data_set.size() - 1); //0λ������Ч����
		}
		static double Summary(const std::vector<double>& data_set) //�������ݼ��ܺ�
		{
			double sum = 0;
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				sum = sum + *vdit;
			}
			return sum;
		}
		static double Summary(const std::vector<double>& data_set, const int power) //�������ݼ���Ԫ��power���ݺ���ܺ�
		{
			double sum = 0;
			std::vector<double> temp = data_set;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				sum = sum + pow(*vdit, power);
			}
			return sum;
		}
		static double CroseSum(const std::vector<double>& x, const std::vector<double>& y) //���������������ݼ���Ӧλ��Ԫ��֮�����ܺͣ���Ϊ�����
		{
			if (Capacity(x) != Capacity(y)) return (0);
			double crssum = 0;
			for (int i = 1; i <= Capacity(x); i++) {
				crssum = crssum + x[i] * y[i];
			}
			return (crssum);
		}
		static double Average(const std::vector<double>& data_set) //�������ݼ���ƽ����
		{
			double sum = Summary(data_set);
			return (sum / (double)(Capacity(data_set)));
		}
		static double CroseAver(const std::vector<double>& data_set1, const std::vector<double>& data_set2) //���������������ݼ��Ľ���ƽ����
		{
			return (CroseSum(data_set1, data_set2) / (double)(Capacity(data_set1)));
		}
		static double SquareAver(const std::vector<double>& data_set) //�������ݼ���Ԫ��ƽ�����ƽ����
		{
			return (Summary(data_set, 2) / (double)(Capacity(data_set)));
		}
		static double StdDev(const std::vector<double>& data_set) //�������ݼ��ı�׼��
		{
			std::vector<double> temp = data_set;
			double aver = Average(temp);
			double dif_sq = 0;
			for (std::vector<double>::iterator vdit = temp.begin() + 1; vdit != temp.end(); vdit++) {
				dif_sq = dif_sq + (aver - *vdit) * (aver - *vdit);
			}
			return (sqrt(dif_sq / (double)(Capacity(data_set))));
		}
		static double AverStdDev(const std::vector<double>& data_set) //�������ݼ���ƽ��ֵ��׼��
		{
			return (StdDev(data_set) / sqrt((double)(Capacity(data_set) - 1)));
		}
		static double MedianNum(const std::vector<double>& data_set) //�������ݼ�����λ��
		{
			std::vector<double> temp = data_set;
			sort(temp.begin() + 1, temp.end()); //���򣬲��Һ���0λ��
			int size = Capacity(temp); //�����������������ż
			if (size % 2 == 0) {
				return ((temp[(size / 2)] + temp[(size / 2) + 1]) / 2.0);
			}
			else return (temp[(size / 2) + 1]);
		}
		static int ModeNum(const std::vector<double>& data_set, std::vector<int>& mode_nums) //�������ݼ�����������mode_nums��װ�ҵ�������������ֵ������һ���������ֵĴ���
		{
			int count = 1;
			struct dat {
				int _num;
				int _count;
			};
			std::vector<dat> stat;
			mode_nums.resize(0);
			std::vector<int> data_set_int(data_set.begin() + 1, data_set.end()); //�����������ݼ��Ǹ�������ת��Ϊ����
			sort(data_set_int.begin(), data_set_int.end()); //����
			for (std::vector<int>::iterator viit = data_set_int.begin() + 1; viit != data_set_int.end(); viit++) {
				if (*(viit) == *(viit - 1)) count++; //�����ǰ����ǰһ������ȣ���ô������count+1
				else {
					stat.push_back({ *(viit - 1), count });
					count = 1; //��������1
				}
				if (viit == data_set_int.end() - 1) stat.push_back({ *(viit - 1), count }); //����ǰ�������һ����������Ҫ��pushһ��
			}
			sort(stat.begin(), stat.end(), [](dat x, dat y) {return x._count > y._count; }); //��dat�е�_count��С��������
			const int max = stat[0]._count;
			for (int i = 0; stat[i]._count == max; i++) mode_nums.push_back(stat[i]._num);
			return stat[0]._count;
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
		static void Print(const std::vector<double>& x) //��/��������
		{
			Process::Print(x);
		}
		static double K(const std::vector<double>& x, const std::vector<double>& y) //�������ֱ�ߵ�б��
		{
			int n = x.size() - 1;
			if (y.size() - 1 != n) return 0;
			else return ((Process::CroseSum(x, y) - (double)n * Process::Average(x) * Process::Average(y)) / (Process::Summary(x, 2) - (double)n * Process::Average(x) * Process::Average(x)));
		}
		static double B(const std::vector<double>& x, const std::vector<double>& y) //�������ֱ�ߵĽؾ�
		{
			return (Process::Average(y) - K(x, y) * Process::Average(x));
		}
		static double R(const std::vector<double>& x, const std::vector<double>& y) //������϶�
		{
			double x_sqaver = Process::SquareAver(x), y_sqaver = Process::SquareAver(y);
			double x_aver = Process::Average(x), y_aver = Process::Average(y);
			return ((Process::CroseAver(x, y) - x_aver * y_aver) / (sqrt((x_sqaver - x_aver * x_aver) * (y_sqaver - y_aver * y_aver))));
		}
	};

};