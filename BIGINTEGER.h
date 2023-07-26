#pragma once
#include <iostream>
#include <vector>
#include <string>

class BigInteger
{
public:
	BigInteger()
	{
		sign = '+';
		_number.clear();
	}
	BigInteger(const long long num)
	{
		_number.clear();
		long long copy = num;
		if (copy < 0) {
			sign = '-';
			copy = -copy;
			while (copy != 0) {
				int temp = copy % 10;
				_number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy > 0) {
			sign = '+';
			while (copy != 0) {
				int temp = copy % 10;
				_number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy == 0) {
			sign = '+';
			_number.push_back(0);
		}
	}
	BigInteger(const std::string& num)
	{
		_number.clear();
		std::string copy = num;
		if ((char)copy[0] == '+') {
			sign = '+';
			copy.erase(copy.begin());
		}
		if ((char)copy[0] == '-') {
			sign = '-';
			copy.erase(copy.begin());
		}
		for (std::string::iterator sit = copy.begin(); sit != copy.end(); sit++) {
			short temp = (short)(*sit - '0');
			_number.push_back(temp);
		}
		std::reverse(_number.begin(), _number.end());
		std::vector<short>::iterator vsit = _number.end() - 1;
		while (*vsit == 0 && vsit != _number.begin()) {
			vsit = vsit - 1;
			_number.erase(vsit + 1);
		}
	}
	BigInteger(const BigInteger& num)
	{
		_number = num._number;
	}
	void Input()
	{
		_number.clear();
		std::string strtemp;
		std::cin >> strtemp;
		if ((char)strtemp[0] == '+') {
			sign = '+';
			strtemp.erase(strtemp.begin());
		}
		if ((char)strtemp[0] == '-') {
			sign = '-';
			strtemp.erase(strtemp.begin());
		}
		for (std::string::iterator sit = strtemp.begin(); sit != strtemp.end(); sit++) {
			short temp = (short)(*sit - '0');
			_number.push_back(temp);
		}
		std::reverse(_number.begin(), _number.end());
		std::vector<short>::iterator vsit = _number.end() - 1;
		while (*vsit == 0 && vsit != _number.begin()) {
			vsit = vsit - 1;
			_number.erase(vsit + 1);
		}
		if (_number.size() == 0) _number.push_back(0);
	}
	void Print()
	{
		if (sign == '+');
		else std::cout << sign;
		for (std::vector<short>::iterator vsit = _number.end() - 1; vsit >= _number.begin() + 1; vsit--) std::cout << *vsit;
		std::cout << *(_number.begin());
	}
	int Digit()
	{
		return this->_number.size();
	}
	BigInteger& operator=(const BigInteger& num)
	{
		this->sign = num.sign;
		this->_number = num._number;
		return (*this);
	}
	BigInteger& operator=(const long long num)
	{
		BigInteger temp(num);
		*this = temp;
		return (*this);
	}
	BigInteger& operator=(const std::string& num)
	{
		BigInteger temp(num);
		*this = temp;
		return (*this);
	}
	BigInteger operator+(const BigInteger& num)
	{
		BigInteger temp = *this, copy = num;
		if (copy.sign == '-') return operator-(-copy); //如果num是负数，那么相当于减去-num
		else { //如果num是正数
			if (temp.sign == '+') { //如果this是正数或0
				temp.sign = '+';
				if (temp.Digit() >= copy.Digit()) {
					while (copy.Digit() < temp.Digit()) copy._number.push_back(0);
				}
				else {
					while (temp.Digit() < copy.Digit()) temp._number.push_back(0);
				}
				for (int i = 0; i < temp.Digit(); i++) temp._number[i] = temp._number[i] + copy._number[i];
				for (int i = 0; i < temp.Digit() - 1; i++) {
					if (temp._number[i] >= 10) {
						temp._number[i + 1] = temp._number[i + 1] + 1;
						temp._number[i] = temp._number[i] - 10;
					}
				}
				int end = temp.Digit() - 1;
				if (temp._number[end] >= 10) {
					temp._number.push_back(temp._number[end] - 10);
					temp._number[end] = 1;
				}
				return temp;
			}
			else { //如果this是负数，那么相当于num - -this
				temp = (copy - (-temp));
				return temp;
			}
		}
	}
	BigInteger& operator-()
	{
		if ((*this).sign == '+') (*this).sign = '-';
		else (*this).sign = '+';
		return (*this);
	}
	BigInteger operator-(const BigInteger& num)
	{
		BigInteger temp = *this, copy = num;
		if (copy.sign == '-') return operator+(-copy); //如果num是负数，那么相当于加上-num
		else { //如果num是正数
			if (temp.sign == '+') { //如果this是正数或0
				if (temp < copy) { //如果this小于num，那么相当于-(num - this)
					temp = (-(copy - temp));
					return temp;
				}
				if (temp.Digit() >= copy.Digit()) {
					while (copy.Digit() < temp.Digit()) copy._number.push_back(0);
				}
				else {
					while (temp.Digit() < copy.Digit()) temp._number.push_back(0);
				}
				for (int i = 0; i < temp.Digit(); i++) temp._number[i] = temp._number[i] - copy._number[i];
				for (int i = 0; i < temp.Digit() - 1; i++) {
					if (temp._number[i] < 0) {
						temp._number[i + 1] = temp._number[i + 1] - 1;
						temp._number[i] = temp._number[i] + 10;
					}
				}
				EraseZero(temp);
				return temp;
			}
			else {
				temp = (-((-temp) + copy)); //如果this是负数，那么相当于-(-this + -num)
				return temp;
			}
		}
	}
	BigInteger operator*(const BigInteger& num)
	{
		BigInteger result;
		if (num.sign == this->sign) {
			result.sign = '+';

		}
		else {
			result.sign = '-';

		}
		return result;
	}
	BigInteger operator/(const BigInteger& num)
	{
		BigInteger result;
		if (num.sign == this->sign) {
			result.sign = '+';

		}
		else {
			result.sign = '-';

		}
		return result;
	}
	bool operator>(const BigInteger& num)
	{
		BigInteger copy = num;
		if (this->sign == '+' && copy.sign == '-') return true;
		if ((*this).Digit() > copy.Digit()) return true;
		if ((*this).Digit() < copy.Digit()) return false;
		if ((*this).Digit() == copy.Digit()) {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] > copy._number[i]) return true;
				if ((*this)._number[i] < copy._number[i]) return false;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
		}
		return false;
	}
	bool operator>=(const BigInteger& num)
	{
		if ((*this) > num || (*this) == num) return true;
		else return false;
	}
	bool operator<(const BigInteger& num)
	{
		BigInteger copy = num;
		if (this->sign == '-' && copy.sign == '+') return true;
		if ((*this).Digit() < copy.Digit()) return true;
		if ((*this).Digit() > copy.Digit()) return false;
		if ((*this).Digit() == copy.Digit()) {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] < copy._number[i]) return true;
				if ((*this)._number[i] > copy._number[i]) return false;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
		}
		return false;
	}
	bool operator<=(const BigInteger& num)
	{
		if ((*this) < num || (*this) == num) return true;
		else return false;
	}
	bool operator==(const BigInteger& num)
	{
		BigInteger copy = num;
		if ((*this).sign != copy.sign || (*this).Digit() != copy.Digit()) return false;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != copy._number[i]) return false;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
			return true;
		}
	}
	bool operator!=(const BigInteger& num)
	{
		BigInteger copy = num;
		if ((*this).sign != copy.sign || (*this).Digit() != copy.Digit()) return true;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != copy._number[i]) return true;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
			return false;
		}
	}

protected:
	char sign = '+';
	std::vector<short> _number;
	void EraseZero(BigInteger& num)
	{
		std::vector<short>::iterator vsit = num._number.end() - 1;
		while (*vsit == 0 && vsit != num._number.begin()) {
			vsit = vsit - 1;
			(*this)._number.erase(vsit + 1);
		}
	}

};