#pragma once
#include <iostream>
#include <vector>
#include <string>

class BigInteger
{
public:
	BigInteger()
	{
		_sign = '+'; //Ĭ��Ϊ������
		_number.clear();
	}
	BigInteger(const long long num)
	{
		_number.clear();
		long long copy = num;
		if (copy < 0) { //���numС��0
			_sign = '-'; //����Ϊ��
			copy = -copy;
			while (copy != 0) {
				int temp = copy % 10; //ÿ��ȡ���һλ��
				_number.push_back(temp); //�������_number
				copy = copy / 10; //ȥ�����һλ��
			}
		}
		else if (copy > 0) { //���num����0
			_sign = '+'; //����Ϊ��
			while (copy != 0) {
				int temp = copy % 10;
				_number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy == 0) { //���num����0
			_sign = '+';
			_number.push_back(0);
		}
	}
	BigInteger(const std::vector<short>& num)
	{
		_sign = '+';
		_number = num;
	}
	BigInteger(const std::string& num)
	{
		_sign = '+';
		_number.clear();
		std::string copy = num;
		if ((char)copy[0] == '+') { //����ַ�����λΪ+
			_sign = '+';
			copy.erase(copy.begin()); //ȥ����λ
		}
		if ((char)copy[0] == '-') { //����ַ�����λΪ-
			_sign = '-';
			copy.erase(copy.begin()); //ȥ����λ
		}
		for (std::string::iterator sit = copy.begin(); sit != copy.end(); sit++) { //������������
			if ((char)*sit == '.') break; //����ַ���������С������ֱ������
			short temp = (short)(*sit - '0');
			_number.push_back(temp);
		}
		std::reverse(_number.begin(), _number.end()); //����洢
		std::vector<short>::iterator vsit = _number.end() - 1;
		while (*vsit == 0 && vsit != _number.begin()) { //ȥ�����������ͷ0
			vsit = vsit - 1;
			_number.erase(vsit + 1);
		}
	}
	BigInteger(const BigInteger& num)
	{
		_number = num._number;
	}
	void Input() //���ַ������캯����ͬ
	{
		_number.clear();
		std::string strtemp;
		std::cin >> strtemp;
		if ((char)strtemp[0] == '+') {
			_sign = '+';
			strtemp.erase(strtemp.begin());
		}
		if ((char)strtemp[0] == '-') {
			_sign = '-';
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
		if (_sign == '+');
		else std::cout << _sign;
		for (std::vector<short>::iterator vsit = temp._number.end() - 1; vsit >= temp._number.begin() + 1; vsit--) std::cout << *vsit;
		std::cout << *(temp._number.begin());
	}
	int Digit() const //���ش�����λ��
	{
		return this->_number.size();
	}
	BigInteger Power(const long long power) const //������power����
	{
		BigInteger result(1), copy = *this;
		if (power > 0) { //��������
			long long pow = power;
			while (pow) { //ʹ�ÿ�����
				if (pow & 1) result = result * copy;
				copy = copy * copy;
				pow = pow / 2;
			}
			return result;
		}
		if (power == 0) return (1); //0����
		if (power < 0) return ((BigInteger)1 / copy).Power(-power); //��������
	}
	BigInteger& operator-() //����ԭ�����෴��
	{
		if ((*this)._sign == '+') (*this)._sign = '-';
		else (*this)._sign = '+';
		return (*this);
	}
	BigInteger& operator=(const BigInteger& num)
	{
		this->_sign = num._sign;
		this->_number = num._number;
		return (*this);
	}
	BigInteger operator+(const BigInteger& num) const
	{
		BigInteger result = *this, copy = num;
		if (copy._sign == '-') return operator-(-copy); //���num�Ǹ�������ô�൱�ڼ�ȥ-num
		else { //���num������
			if (result._sign == '+') { //���this��������0
				result._sign = '+';
				if (result.Digit() >= copy.Digit()) while (copy.Digit() < result.Digit()) copy._number.push_back(0); //�������λ�����ڵ���num��numǰ��0�Զ���
				else while (result.Digit() < copy.Digit()) result._number.push_back(0); //�������λ��С��num������ǰ��0�Զ���
				for (int i = 0; i < result.Digit(); i++) result._number[i] = result._number[i] + copy._number[i]; //����Ӧλ��ֱ�����
				for (int i = 0; i < result.Digit() - 1; i++) { //��λ
					if (result._number[i] >= 10) {
						result._number[i + 1] = result._number[i + 1] + 1;
						result._number[i] = result._number[i] - 10;
					}
				}
				int end = result.Digit() - 1; //ȡ���λ
				if (result._number[end] >= 10) { //������λ����������10
					result._number.push_back(1); //�����λ���һλ��Ϊ1
					result._number[end] = result._number[end] - 10; //��ǰλȡ��λ��
				}
				return result;
			}
			else return (copy - (-result)); //���this�Ǹ�������ô�൱��num - -this
		}
	}
	BigInteger operator-(const BigInteger& num) const
	{
		BigInteger result = *this, copy = num;
		if (copy._sign == '-') return operator+(-copy); //���num�Ǹ�������ô�൱�ڼ���-num
		else { //���num������
			if (result._sign == '+') { //���this��������0
				if (result < copy) return (-(copy - result)); //���thisС��num����ô�൱��-(num - this)
				if (result.Digit() >= copy.Digit()) while (copy.Digit() < result.Digit()) copy._number.push_back(0); //�������λ�����ڵ���num��numǰ��0�Զ���
				else while (result.Digit() < copy.Digit()) result._number.push_back(0); //�������λ��С��num������ǰ��0�Զ���
				for (int i = 0; i < result.Digit(); i++) result._number[i] = result._number[i] - copy._number[i]; //����Ӧλ��ֱ�����
				for (int i = 0; i < result.Digit() - 1; i++) { //��λ
					if (result._number[i] < 0) {
						result._number[i + 1] = result._number[i + 1] - 1;
						result._number[i] = result._number[i] + 10;
					}
				}
				EraseZero(result); //Ĩȥͷ0
				return result;
			}
			else return (-((-result) + copy)); //���this�Ǹ�������ô�൱��-(-this + -num)
		}
	}
	BigInteger operator*(const BigInteger& num) const
	{
		BigInteger result, copy = *this;
		result._number.resize(this->Digit() + num.Digit(), 0); //resultλ�����������λ��֮��
		if (num._sign == this->_sign) result._sign = '+'; //�������������ͬresult����Ϊ+
		else result._sign = '-'; //���������������result����Ϊ-
		for (int i = 0, numdigit = num.Digit(); i < numdigit; i++)
			for (int j = 0, thisdigit = this->Digit(); j < thisdigit; j++) result._number[i + j] = result._number[i + j] + this->_number[j] * num._number[i];
		for (int i = 0, redigit = result.Digit(); i < redigit; i++) { //��λ
			BigInteger temp(result._number[i]);
			for (int j = 0; j < temp.Digit(); j++) {
				if (j == 0) result._number[i] = temp._number[0];
				else result._number[i + j] = result._number[i + j] + temp._number[j];
			}
		}
		EraseZero(result);
		return result;
	}
	BigInteger operator/(const BigInteger& num) const
	{
		BigInteger result, remain;
		result = Division(*this, num, remain);
		if (num._sign == this->_sign) result._sign = '+';
		else result._sign = '-';
		return result;
	}
	BigInteger operator%(const BigInteger& num) const
	{
		BigInteger remain;
		Division(*this, num, remain);
		return remain;
	}
	bool operator>(const BigInteger& num) const
	{
		BigInteger copy = num;
		if (this->_sign == '+' && copy._sign == '-') return true;
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
	bool operator>=(const BigInteger& num) const
	{
		if (this->_sign == '+' && num._sign == '-') return true;
		if ((*this).Digit() > num.Digit()) return true;
		if ((*this).Digit() < num.Digit()) return false;
		if ((*this).Digit() == num.Digit()) {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] > num._number[i]) return true;
				if ((*this)._number[i] < num._number[i]) return false;
				if ((*this)._number[i] == num._number[i]) continue;
			}
		}
		return true;
	}
	bool operator<(const BigInteger& num) const
	{
		BigInteger copy = num;
		if (this->_sign == '-' && copy._sign == '+') return true;
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
	bool operator<=(const BigInteger& num) const
	{
		if (this->_sign == '-' && num._sign == '+') return true;
		if ((*this).Digit() < num.Digit()) return true;
		if ((*this).Digit() > num.Digit()) return false;
		if ((*this).Digit() == num.Digit()) {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] < num._number[i]) return true;
				if ((*this)._number[i] > num._number[i]) return false;
				if ((*this)._number[i] == num._number[i]) continue;
			}
		}
		return true;
	}
	bool operator==(const BigInteger& num) const
	{
		BigInteger copy = num;
		if ((*this)._sign != copy._sign || (*this).Digit() != copy.Digit()) return false;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != copy._number[i]) return false;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
			return true;
		}
	}
	bool operator!=(const BigInteger& num) const
	{
		BigInteger copy = num;
		if ((*this)._sign != copy._sign || (*this).Digit() != copy.Digit()) return true;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != copy._number[i]) return true;
				if ((*this)._number[i] == copy._number[i]) continue;
			}
			return false;
		}
	}
	explicit operator long long()
	{
		int min = (15 < this->Digit() - 1) ? 15 : this->Digit() - 1;
		long long result = 0;
		for (int i = 0; i <= min; i++) result = result + this->_number[i] * pow(10, i);
		if (this->_sign == '-') result = -result;
		return result;
	}

	static BigInteger Division(const BigInteger& dividend, const BigInteger& divisor, BigInteger& remainder) //dividend / divisor = quotient...remainder
	{
		if (divisor == (long long)0) exit(0);
		if (dividend == (long long)0 || dividend < divisor) {
			remainder = dividend;
			return (long long)0;
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
	char _sign = '+';
	std::vector<short> _number;
	static void AddZero(BigInteger& num) //��10ʱ�൱��ֱ��ĩλ��0
	{
		num._number.insert(num._number.begin(), 0);
	}
	static void AddZero(BigInteger& num, int power) //��10^powerʱ�൱��ֱ��ĩλ��power��0
	{
		for (int i = 1; i <= power; i++) num._number.insert(num._number.begin(), 0);
	}
	static void EraseZero(BigInteger& num) //Ĩ��ͷ0
	{
		std::vector<short>::iterator vsit = num._number.end() - 1;
		while (*vsit == 0 && vsit != num._number.begin()) {
			vsit = vsit - 1;
			num._number.erase(vsit + 1);
		}
	}

};