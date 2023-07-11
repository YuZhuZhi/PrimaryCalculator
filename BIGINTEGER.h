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
		number.clear();
	}
	BigInteger(const long long num)
	{
		number.clear();
		long long copy = num;
		if (copy < 0) {
			sign = '-';
			copy = -copy;
			while (copy != 0) {
				int temp = copy % 10;
				number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy > 0) {
			sign = '+';
			while (copy != 0) {
				int temp = copy % 10;
				number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy == 0) {
			sign = '+';
			number.push_back(0);
		}
	}
	BigInteger(const std::string& num)
	{
		number.clear();
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
			number.push_back(temp);
		}
		std::reverse(number.begin(), number.end());
		std::vector<short>::iterator vsit = number.end() - 1;
		while (*vsit == 0 && vsit != number.begin()) {
			vsit = vsit - 1;
			number.erase(vsit + 1);
		}
	}
	void Input()
	{
		number.clear();
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
			number.push_back(temp);
		}
		std::reverse(number.begin(), number.end());
		std::vector<short>::iterator vsit = number.end() - 1;
		while (*vsit == 0 && vsit != number.begin()) {
			vsit = vsit - 1;
			number.erase(vsit + 1);
		}
		if (number.size() == 0) number.push_back(0);
	}
	void Print()
	{
		if (sign == '+');
		else std::cout << sign;
		for (std::vector<short>::iterator vsit = number.end() - 1; vsit >= number.begin() + 1; vsit--) std::cout << *vsit;
		std::cout << *(number.begin());
	}
	inline int Digit()
	{
		return this->number.size();
	}
	BigInteger& operator=(const BigInteger& num)
	{
		this->sign = num.sign;
		this->number = num.number;
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
	BigInteger& operator+(const BigInteger& num)
	{
		BigInteger copy = num;
		if (copy.sign == '-') return operator-(-copy); //如果num是负数，那么相当于减去-num
		else { //如果num是正数
			if (this->sign == '+') { //如果this是正数或0
				this->sign = '+';
				if ((*this).Digit() >= copy.Digit()) {
					while (copy.Digit() < (*this).Digit()) copy.number.push_back(0);
				}
				else {
					while ((*this).Digit() < copy.Digit()) (*this).number.push_back(0);
				}
				for (int i = 0; i < (*this).Digit(); i++) this->number[i] = this->number[i] + copy.number[i];
				for (int i = 0; i < (*this).Digit() - 1; i++) {
					if (this->number[i] >= 10) {
						this->number[i + 1] = this->number[i + 1] + 1;
						this->number[i] = this->number[i] - 10;
					}
				}
				int end = (*this).Digit() - 1;
				if (this->number[end] >= 10) {
					this->number.push_back(this->number[end] - 10);
					this->number[end] = 1;
				}
				return (*this);
			}
			else { //如果this是负数，那么相当于num - -this
				(*this) = (copy - (-(*this)));
				return (*this);
			}
		}
	}
	inline BigInteger& operator-()
	{
		if ((*this).sign == '+') (*this).sign = '-';
		else (*this).sign = '+';
		return (*this);
	}
	BigInteger& operator-(const BigInteger& num)
	{
		BigInteger copy = num;
		if (copy.sign == '-') return operator+(-copy); //如果num是负数，那么相当于加上-num
		else { //如果num是正数
			if ((*this).sign == '+') { //如果this是正数或0
				if ((*this) < copy) { //如果this小于num，那么相当于-(num - this)
					(*this) = (-(copy - *this));
					return (*this);
				}
				if ((*this).Digit() >= copy.Digit()) {
					while (copy.Digit() < (*this).Digit()) copy.number.push_back(0);
				}
				else {
					while ((*this).Digit() < copy.Digit()) (*this).number.push_back(0);
				}
				for (int i = 0; i < (*this).Digit(); i++) this->number[i] = this->number[i] - copy.number[i];
				for (int i = 0; i < (*this).Digit() - 1; i++) {
					if (this->number[i] < 0) {
						this->number[i + 1] = this->number[i + 1] - 1;
						this->number[i] = this->number[i] + 10;
					}
				}
				EraseZero(*this);
				return (*this);
			}
			else {
				(*this) = (-((-(*this)) + copy)); //如果this是负数，那么相当于-(-this + -num)
				return (*this);
			}
		}
	}
	BigInteger& operator*(const BigInteger& num)
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
	BigInteger& operator/(const BigInteger& num)
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
				if ((*this).number[i] > copy.number[i]) return true;
				if ((*this).number[i] < copy.number[i]) return false;
				if ((*this).number[i] == copy.number[i]) continue;
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
				if ((*this).number[i] < copy.number[i]) return true;
				if ((*this).number[i] > copy.number[i]) return false;
				if ((*this).number[i] == copy.number[i]) continue;
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
				if ((*this).number[i] != copy.number[i]) return false;
				if ((*this).number[i] == copy.number[i]) continue;
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
				if ((*this).number[i] != copy.number[i]) return true;
				if ((*this).number[i] == copy.number[i]) continue;
			}
			return false;
		}
	}

protected:
	char sign = '+';
	std::vector<short> number;
	void EraseZero(BigInteger& num)
	{
		std::vector<short>::iterator vsit = num.number.end() - 1;
		while (*vsit == 0 && vsit != num.number.begin()) {
			vsit = vsit - 1;
			(*this).number.erase(vsit + 1);
		}
	}
};