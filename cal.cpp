#define _CRT_SECURE_NO_WARNINGS
#include <vector>
#include <Windows.h>
#include "DISPLAY.h"
#include "LINEAR.h"
using namespace std;

char YN = 'N';
int choicefir = 0, choicesec = 0, choicethr = 0;
vector<double> data_set, x_set, y_set, right_value;
Linear::Vector<double> x, y;
Linear::Matrix<double> matrix, matrix1, matrix2, matrix3, coeff;

int main()
{
	system("color F0");
	Display::Guide();
	while (true) {
		Display::ChooseFir();
	}
}