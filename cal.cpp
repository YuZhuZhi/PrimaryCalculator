#define _CRT_SECURE_NO_WARNINGS
#include <cstdio>
#include <iostream>
#include <vector>
#include <Windows.h>
#include "DATA.h"
#include "LINEAR.h"
#include "DISPLAY.h"
using namespace std;

char YN = 'N';
int choicefir = 0, choicesec = 0, choicethr = 0;
vector<double> data_set, x, y, right_value;
vector<vector<double> > matrix, matrix1, matrix2, matrix3, coeff;

int main()
{
	system("color F0");
	Display::Guide();
	while (true) {
		Display::ChooseFir();
	}
}