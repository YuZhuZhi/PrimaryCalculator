#pragma once
#include <iostream>
#include <vector>
#include <string>

class BigInteger
{
public:
	BigInteger()
	{
		_sign = '+'; //默认为正符号
		_number.clear();
	}
	BigInteger(const long long num)
	{
		_number.clear();
		long long copy = num;
		if (copy < 0) { //如果num小于0
			_sign = '-'; //符号为负
			copy = -copy;
			while (copy != 0) {
				int temp = copy % 10; //每次取最后一位数
				_number.push_back(temp); //倒序插入_number
				copy = copy / 10; //去掉最后一位数
			}
		}
		else if (copy > 0) { //如果num大于0
			_sign = '+'; //符号为正
			while (copy != 0) {
				int temp = copy % 10;
				_number.push_back(temp);
				copy = copy / 10;
			}
		}
		else if (copy == 0) { //如果num等于0
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
		if ((char)copy[0] == '+') { //如果字符串首位为+
			_sign = '+';
			copy.erase(copy.begin()); //去掉首位
		}
		if ((char)copy[0] == '-') { //如果字符串首位为-
			_sign = '-';
			copy.erase(copy.begin()); //去掉首位
		}
		for (std::string::iterator sit = copy.begin(); sit != copy.end(); sit++) { //复制整数部分
			if ((char)*sit == '.') break; //如果字符串数字有小数部分直接舍弃
			int temp = (int)(*sit - '0');
			_number.push_back(temp);
		}
		std::reverse(_number.begin(), _number.end()); //倒序存储
		auto vsit = _number.end() - 1;
		while (*vsit == 0 && vsit != _number.begin()) { //去除恶意输入的头0
			vsit = vsit - 1;
			_number.erase(vsit + 1);
		}
	}
	BigInteger(const BigInteger& num)
	{
		_number = num._number;
	}
	void Input() //与字符串构造函数相同
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
		auto vsit = _number.end() - 1;
		while (*vsit == 0 && vsit != _number.begin()) {
			vsit = vsit - 1;
			_number.erase(vsit + 1);
		}
		if (_number.size() == 0) _number.push_back(0);
	}
	void Print() const
	{
		if (_sign == '-') std::cout << '-';
		for (auto vsit = this->_number.end() - 1; vsit >= this->_number.begin() + 1; vsit--) std::cout << *vsit;
		std::cout << *(this->_number.begin());
	}
	int Digit() const //返回大数的位数
	{
		return (this->_number.size());
	}
	BigInteger Power(const long long power) const //大数的power次幂
	{
		BigInteger result(1), copy = *this;
		if (power > 0) { //正数次幂
			long long pow = power;
			while (pow) { //使用快速幂
				if (pow & 1) result = result * copy;
				copy = copy * copy;
				pow = pow / 2;
			}
			return result;
		}
		if (power == 0) return (1); //0次幂
		if (power < 0) return ((BigInteger)1 / copy).Power(-power); //负数次幂
	}
	BigInteger& operator-() //返回原数的相反数
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
		BigInteger result(*this), copy(num);
		if (copy._sign == '-') return operator-(-copy); //如果num是负数，那么相当于减去-num
		else { //如果num是正数
			if (result._sign == '+') { //如果this是正数或0
				result._sign = '+';
				if (result.Digit() >= copy.Digit()) while (copy.Digit() < result.Digit()) copy._number.push_back(0); //如果此数位数大于等于num，num前补0以对齐
				else while (result.Digit() < copy.Digit()) result._number.push_back(0); //如果此数位数小于num，此数前补0以对齐
				for (int i = 0; i < result.Digit(); i++) result._number[i] = result._number[i] + copy._number[i]; //按相应位数直接相加
				for (int i = 0; i < result.Digit() - 1; i++) { //进位
					if (result._number[i] >= 10) {
						result._number[i + 1] = result._number[i + 1] + 1;
						result._number[i] = result._number[i] - 10;
					}
				}
				int end = result.Digit() - 1; //取最高位
				if (result._number[end] >= 10) { //如果最高位所存数大于10
					result._number.push_back(1); //向最高位添加一位数为1
					result._number[end] = result._number[end] - 10; //当前位取个位数
				}
				return result;
			}
			else return (copy - (-result)); //如果this是负数，那么相当于num - -this
		}
	}
	BigInteger operator-(const BigInteger& num) const
	{
		BigInteger result(*this), copy(num);
		if (copy._sign == '-') return operator+(-copy); //如果num是负数，那么相当于加上-num
		else { //如果num是正数
			if (result._sign == '+') { //如果this是正数或0
				if (result < copy) return (-(copy - result)); //如果this小于num，那么相当于-(num - this)
				if (result.Digit() >= copy.Digit()) while (copy.Digit() < result.Digit()) copy._number.push_back(0); //如果此数位数大于等于num，num前补0以对齐
				else while (result.Digit() < copy.Digit()) result._number.push_back(0); //如果此数位数小于num，此数前补0以对齐
				for (int i = 0; i < result.Digit(); i++) result._number[i] = result._number[i] - copy._number[i]; //按相应位数直接相减
				for (int i = 0; i < result.Digit() - 1; i++) { //借位
					if (result._number[i] < 0) {
						result._number[i + 1] = result._number[i + 1] - 1;
						result._number[i] = result._number[i] + 10;
					}
				}
				result.EraseZero(); //抹去头0
				return result;
			}
			else return (-((-result) + copy)); //如果this是负数，那么相当于-(-this + -num)
		}
	}
	BigInteger operator*(const BigInteger& num) const
	{
		BigInteger result;
		result._number.resize(this->Digit() + num.Digit(), 0); //result位数最多是两数位数之和
		if (num._sign == this->_sign) result._sign = '+'; //如果两数符号相同result符号为+
		else result._sign = '-'; //如果两数符号相异result符号为-
		for (int i = 0, numdigit = num.Digit(); i < numdigit; i++)
			for (int j = 0, thisdigit = this->Digit(); j < thisdigit; j++) result._number[i + j] = result._number[i + j] + this->_number[j] * num._number[i];
		for (int i = 0, redigit = result.Digit(); i < redigit; i++) { //进位
			BigInteger temp(result._number[i]);
			for (int j = 0; j < temp.Digit(); j++) {
				if (j != 0) result._number[i + j] = result._number[i + j] + temp._number[j];
				else result._number[i] = temp._number[0];
			}
		}
		result.EraseZero();
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
		if ((*this)._sign != num._sign || (*this).Digit() != num.Digit()) return false;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != num._number[i]) return false;
				if ((*this)._number[i] == num._number[i]) continue;
			}
			return true;
		}
	}
	bool operator!=(const BigInteger& num) const
	{
		if ((*this)._sign != num._sign || (*this).Digit() != num.Digit()) return true;
		else {
			int end = (*this).Digit() - 1;
			for (int i = end; i >= 0; i--) {
				if ((*this)._number[i] != num._number[i]) return true;
				if ((*this)._number[i] == num._number[i]) continue;
			}
			return false;
		}
	}
	explicit operator long long() const
	{
		int min = (19 < this->Digit() - 1) ? 19 : this->Digit() - 1;
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
			temp.EraseZero();
			int count = 0;
			while (temp >= divisor) {
				temp = temp - divisor;
				count++;
			}
			quotient._number.insert(quotient._number.begin(), count);
		}
		remainder = temp;
		quotient.EraseZero();
		return quotient;
	}

protected:
	char _sign = '+';
	std::vector<short> _number;
	void AddZero() //乘10时相当于直接末位加0
	{
		this->_number.insert(this->_number.begin(), 0);
	}
	void AddZero(int power) //乘10^power时相当于直接末位加power个0
	{
		for (int i = 1; i <= power; i++) this->_number.insert(this->_number.begin(), 0);
	}
	void EraseZero() //抹除头0
	{
		auto vsit = this->_number.end() - 1;
		while (*vsit == 0 && vsit != this->_number.begin()) {
			vsit = vsit - 1;
			this->_number.erase(vsit + 1);
		}
	}

};