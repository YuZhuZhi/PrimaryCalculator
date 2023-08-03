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
	BigInteger(const std::vector<short> num)
	{
		sign = '+';
		_number = num;
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
	void Print() const
	{
		BigInteger temp = *this;
		if (sign == '+');
		else std::cout << sign;
		for (std::vector<short>::iterator vsit = temp._number.end() - 1; vsit >= temp._number.begin() + 1; vsit--) std::cout << *vsit;
		std::cout << *(temp._number.begin());
	}
	int Digit() const
	{
		return this->_number.size();
	}
	BigInteger Power(const long long power) const
	{
		BigInteger result(1), copy = *this;
		if (power > 0) {
			long long pow = power;
			while (pow) {
				if (pow & 1) result = result * copy;
				copy = copy * copy;
				pow = pow / 2;
			}
			return result;
		}
		if (power == 0) return (1);
		if (power < 0) return ((BigInteger)1 / copy).Power(-power);
	}
	BigInteger& operator-()
	{
		if ((*this).sign == '+') (*this).sign = '-';
		else (*this).sign = '+';
		return (*this);
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
	BigInteger operator+(const BigInteger& num) const
	{
		BigInteger temp = *this, copy = num;
		if (copy.sign == '-') return operator-(-copy); //���num�Ǹ�������ô�൱�ڼ�ȥ-num
		else { //���num������
			if (temp.sign == '+') { //���this��������0
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
					temp._number.push_back(1);
					temp._number[end] = temp._number[end] - 10;
				}
				return temp;
			}
			else { //���this�Ǹ�������ô�൱��num - -this
				temp = (copy - (-temp));
				return temp;
			}
		}
	}
	BigInteger operator+(const long long num) const
	{
		return (*this + (BigInteger)num);
	}
	BigInteger operator+(const std::string& num) const
	{
		return (*this + (BigInteger)num);
	}
	BigInteger operator-(const BigInteger& num) const
	{
		BigInteger temp = *this, copy = num;
		if (copy.sign == '-') return operator+(-copy); //���num�Ǹ�������ô�൱�ڼ���-num
		else { //���num������
			if (temp.sign == '+') { //���this��������0
				if (temp < copy) { //���thisС��num����ô�൱��-(num - this)
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
				temp = (-((-temp) + copy)); //���this�Ǹ�������ô�൱��-(-this + -num)
				return temp;
			}
		}
	}
	BigInteger operator-(const long long num) const
	{
		return (*this - (BigInteger)num);
	}
	BigInteger operator-(const std::string& num) const
	{
		return (*this - (BigInteger)num);
	}
	BigInteger operator*(const BigInteger& num) const
	{
		BigInteger result, copy = *this;
		result._number.resize(this->Digit() + num.Digit(), 0);
		if (num.sign == this->sign) {
			result.sign = '+';
			for (int i = 0; i < num.Digit(); i++) {
				for (int j = 0; j < this->Digit(); j++) {
					result._number[i + j] = result._number[i + j] + this->_number[j] * num._number[i];
				}
			}
			for (int i = 0; i < result.Digit(); i++) { //��λ
				BigInteger temp(result._number[i]);
				for (int j = 0; j < temp.Digit(); j++) {
					if (j == 0) result._number[i] = temp._number[0];
					else result._number[i + j] = result._number[i + j] + temp._number[j];
				}
			}
			EraseZero(result);
			return result;
		}
		else {
			result.sign = '-';
			for (int i = 0; i < num.Digit(); i++) {
				for (int j = 0; j < this->Digit(); j++) {
					result._number[i + j] = result._number[i + j] + this->_number[j] * num._number[i];
				}
			}
			for (int i = 0; i < result.Digit(); i++) { //��λ
				BigInteger temp(result._number[i]);
				for (int j = 0; j < temp.Digit(); j++) {
					if (j == 0) result._number[i] = temp._number[0];
					else result._number[i + j] = result._number[i + j] + temp._number[j];
				}
			}
			EraseZero(result);
			return result;
		}
	}
	BigInteger operator*(const long long num) const
	{
		return (*this * (BigInteger)num);
	}
	BigInteger operator*(const std::string& num) const
	{
		return (*this * (BigInteger)num);
	}
	BigInteger operator/(const BigInteger& num) const
	{
		BigInteger result, remain;
		if (num.sign == this->sign) {
			result = Division(*this, num, remain);
			result.sign = '+';
			return result;
		}
		else {
			result = Division(*this, num, remain);
			result.sign = '-';
			return result;
		}
	}
	BigInteger operator/(const long long num) const
	{
		return (*this / (BigInteger)num);
	}
	BigInteger operator/(const std::string& num) const
	{
		return (*this / (BigInteger)num);
	}
	BigInteger operator%(const BigInteger& num) const
	{
		BigInteger result, remain;
		if (num.sign == this->sign) {
			result.sign = '+';
			result = Division(*this, num, remain);
			return remain;
		}
		else {
			result.sign = '-';
			result = Division(*this, num, remain);
			return remain;
		}
	}
	BigInteger operator%(const long long num) const
	{
		return (*this % (BigInteger)num);
	}
	BigInteger operator%(const std::string& num) const
	{
		return (*this % (BigInteger)num);
	}
	bool operator>(const BigInteger& num) const
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
	bool operator>(const long long num) const
	{
		return (*this > (BigInteger)num);
	}
	bool operator>(const std::string& num) const
	{
		return (*this > (BigInteger)num);
	}
	bool operator>=(const BigInteger& num) const
	{
		if ((*this) > num || (*this) == num) return true;
		else return false;
	}
	bool operator>=(const long long num) const
	{
		return (*this >= (BigInteger)num);
	}
	bool operator>=(const std::string& num) const
	{
		return (*this >= (BigInteger)num);
	}
	bool operator<(const BigInteger& num) const
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
	bool operator<(const long long num) const
	{
		return (*this < (BigInteger)num);
	}
	bool operator<(const std::string& num) const
	{
		return (*this < (BigInteger)num);
	}
	bool operator<=(const BigInteger& num) const
	{
		if ((*this) < num || (*this) == num) return true;
		else return false;
	}
	bool operator<=(const long long num) const
	{
		return (*this <= (BigInteger)num);
	}
	bool operator<=(const std::string& num) const
	{
		return (*this <= (BigInteger)num);
	}
	bool operator==(const BigInteger& num) const
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
	bool operator==(const long long num) const
	{
		return (*this == (BigInteger)num);
	}
	bool operator==(const std::string& num) const
	{
		return (*this == (BigInteger)num);
	}
	bool operator!=(const BigInteger& num) const
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
	bool operator!=(const long long num) const
	{
		return (*this != (BigInteger)num);
	}
	bool operator!=(const std::string& num) const
	{
		return (*this != (BigInteger)num);
	}
	operator double()
	{
		int min = (15 < this->Digit() - 1) ? 15 : this->Digit() - 1;
		double result = 0;
		for (int i = 0; i <= min; i++) result = result + this->_number[i] * pow(10, i);
		if (this->sign == '-') result = -result;
		return result;
	}

	static BigInteger Division(const BigInteger& dividend, const BigInteger& divisor, BigInteger& remainder) //dividend / divisor = quotient...remainder
	{
		if (divisor == 0) exit(0);
		if (dividend == 0 || dividend < divisor) {
			remainder = dividend;
			return 0;
		}
		BigInteger quotient, temp;
		for (int i = 0; i < dividend.Digit(); i++) {
			temp._number.insert(temp._number.begin(), dividend._number[dividend.Digit() - 1 - i]);
			EraseZero(temp);
			int count = 0;
			while (temp >= divisor) {
				temp = temp - divisor;
				count++;
			}
			quotient._number.insert(quotient._number.begin(), count);
		}
		remainder = temp;
		EraseZero(quotient);
		return quotient;
	}

protected:
	char sign = '+';
	std::vector<short> _number;
	static void AddZero(BigInteger& num)
	{
		num._number.insert(num._number.begin(), 0);
	}
	static void AddZero(BigInteger& num, int power)
	{
		for (int i = 1; i <= power; i++) num._number.insert(num._number.begin(), 0);
	}
	static void EraseZero(BigInteger& num)
	{
		std::vector<short>::iterator vsit = num._number.end() - 1;
		while (*vsit == 0 && vsit != num._number.begin()) {
			vsit = vsit - 1;
			num._number.erase(vsit + 1);
		}
	}

};