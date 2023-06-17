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
		std::cout << "��ӭʹ����ѧ������ver0.56(TEST)��(Programmed by ����֮)" << std::endl;
		std::cout << "(����֮ ֪����ҳ��https://www.zhihu.com/people/yuzhuzhi)" << std::endl;
		std::cout << "(����������Դ�ڣ�https://github.com/YuZhuZhi/PrimaryCalculator)" << std::endl;
		std::cout << "ʹ�÷�����������ʾ�������ֻ�����ĸ��س���ѡ�����(��ĸ�ɲ����ִ�Сд)��" << std::endl;
		std::cout << "��������ʱ�����ÿո�������ÿ����һ�μ���һ�λس����ύ��" << std::endl;
		std::cout << "Ҳ���Դ�Excel�и������ݺ��ύ��" << std::endl;
		std::cout << "�������뷽�����ύ�����һ�����ݺ󣬰�Ctrl+Z����\"^Z\"���ύ���ɡ�ʾ����" << std::endl;
		std::cout << "1 2 3 4 5" << std::endl << "^Z" << std::endl;
		std::cout << "�̳��ϴ�������ָ����ʹ���ϴ����������ݣ��ٴ�����һ�ֻ���ͬһ�ֲ����С�" << std::endl;
		std::cout << "1e5��ʾ1��10��5�η����Դ����ơ�" << std::endl;
		std::cout << "ף��ʹ����졣" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
	}
	static void ChooseYN()
	{
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cout << "�Ƿ�̳���һ����������ݣ�" << std::endl;
		std::cout << "Y\t��" << std::endl;
		std::cout << "N\t��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����" << std::endl;
		std::cout << "1\t�������ݴ���(6)" << std::endl;
		std::cout << "2\t��ѧ���������ʽ(5)(��δ����)" << std::endl;
		std::cout << "3\t���Դ���(3)" << std::endl << std::endl;
		std::cout << "4\t�������ʾ����" << std::endl;
		std::cout << "--------------------------------------------------------" << std::endl;
		std::cin >> choicefir;
		if (std::cin.fail()) {
			std::cin.clear();
			std::cin.ignore();
			std::cout << "��������С���?" << std::endl;
			std::cout << "������ѡ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����" << std::endl;
		std::cout << "�������ݴ���" << std::endl;
		std::cout << "1\t���" << std::endl;
		std::cout << "2\t���" << std::endl;
		std::cout << "3\t���׼��(��) / ʵ���׼��(��)" << std::endl;
		std::cout << "4\t��ƽ��ֵ��׼��(��) / ƽ��ֵ��ʵ���׼��(��)" << std::endl;
		std::cout << "5\t��ֱ�Ӳ������Ĳ�ȷ����" << std::endl;
		std::cout << "6\t�������������" << std::endl << std::endl;
		std::cout << "7...\t�ص���һ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����" << std::endl;
		std::cout << "��ѧ���������ʽ��" << std::endl;
		std::cout << "1\t�׳�" << std::endl;
		std::cout << "2\t˫�׳�" << std::endl;
		std::cout << "3\t36���ڽ���ת��" << std::endl;
		std::cout << "4\t��չ��������ʽ֮��" << std::endl;
		std::cout << "5\t����ȫ�ֽ�Ķ���ʽ��չ��" << std::endl << std::endl;
		std::cout << "6...\t�ص���һ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����(��������ʱ��������Ctrl+Z��ȷ��)" << std::endl;
		std::cout << "���Դ���:" << std::endl;
		std::cout << "1\t�ⷽ��(��)(4)" << std::endl;
		std::cout << "2\t����(4)" << std::endl;
		std::cout << "3\t����ʽ�����(10)" << std::endl << std::endl;
		std::cout << "4...\t�ص���һ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����(����ʱ��������Ctrl+Z��ȷ��)" << std::endl;
		std::cout << "�ⷽ��(��):" << std::endl;
		std::cout << "1\t��һԪ���η���" << std::endl;
		std::cout << "2\t��һԪ���η���(��δ����)" << std::endl;
		std::cout << "3\t��һԪn�η���(��δ����)" << std::endl;
		std::cout << "4\t��nԪ������" << std::endl << std::endl;
		std::cout << "5...\t�ص���һ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����(����ʱ��������Ctrl+Z��ȷ��)" << std::endl;
		std::cout << "������" << std::endl;
		std::cout << "1\t�������ĳ���" << std::endl;
		std::cout << "2\t��������ı�����" << std::endl;
		std::cout << "3\t���������������" << std::endl;
		std::cout << "4\t��������ļн�(���ػ���)" << std::endl << std::endl;
		std::cout << "5...\t�ص���һ��" << std::endl;
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
		std::cout << "ѡ������Ҫ�Ĳ�����(����ʱ��������Ctrl+Z��ȷ��)" << std::endl;
		std::cout << "(����ʽ/����ɻ���̳У������ɼ̳е����󣻷�֮��Ȼ)" << std::endl;
		std::cout << "����ʽ�����:" << std::endl;
		std::cout << "1\t������ʽ��ֵ" << std::endl;
		std::cout << "2\t�������֮��" << std::endl;
		std::cout << "3\t�������֮��" << std::endl;
		std::cout << "4\t�������֮��" << std::endl;
		std::cout << "5\t�������н����ξ���(��˹��Ԫ��)" << std::endl;
		std::cout << "6\t������������ξ���" << std::endl;
		std::cout << "7\t�����n����" << std::endl;
		std::cout << "8\t����İ�����" << std::endl;
		std::cout << "9\t���������" << std::endl;
		std::cout << "10\t���������ֵ(����)" << std::endl;
		std::cout << "11\t����ʩ����������(��δ����)" << std::endl << std::endl;
		std::cout << "12...\t�ص���һ��" << std::endl;
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