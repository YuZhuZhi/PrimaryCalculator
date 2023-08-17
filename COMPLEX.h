#pragma once
#include <iostream>
#include <cmath>

const double e = 2.71828182846;
const double PI = 3.14159265359;

class Complex
{
public:
	Complex()
	{
		_re = 0;
		_im = 0;
	}
	Complex(const double real)
	{
		_re = real;
		_im = 0;
	}
	Complex(const double real, const double imaginary)
	{
		_re = real;
		_im = imaginary;
	}
	Complex(const double r, const double theta, const bool)
	{
		_re = r * cos(AngleReduct(theta));
		_im = r * sin(AngleReduct(theta));
	}
	Complex(const Complex& complex)
	{
		_re = complex._re;
		_im = complex._im;
	}
	void Input()
	{
		std::cin >> _re;
		std::cin >> _im;
	}
	void Print() const
	{
		if (_re != 0 && _im != 0) {
			if (_im > 0) {
				if (_im != 1) {
					std::cout << _re << " + " << _im << "i" << std::endl;
					return;
				}
				if (_im == 1) {
					std::cout << _re << " + " << "i" << std::endl;
					return;
				}
			}
			if (_im < 0) {
				if (_im != -1) {
					std::cout << _re << " - " << -_im << "i" << std::endl;
					return;
				}
				if (_im == -1) {
					std::cout << _re << " - " << "i" << std::endl;
					return;
				}
			}
		}
		if (_re != 0 && _im == 0) {
			std::cout << _re << std::endl;
			return;
		}
		if (_re == 0 && _im != 0) {
			if (_im != 1 && _im != -1) {
				std::cout << _im << "i" << std::endl;
				return;
			}
			if (_im == 1) {
				std::cout << "i" << std::endl;
				return;
			}
			if (_im == -1) {
				std::cout << "-i" << std::endl;
				return;
			}
		}
		if (_re == 0 && _im == 0) {
			std::cout << 0 << std::endl;
			return;
		}
	}
	void PrintPolar() const
	{
		std::cout << this->Modulus() << "e^{" << Arg(*this) << "i}";
	}
	Complex Unit()
	{
		Complex temp(0, 1);
		return temp;
	}
	Complex Unit(const double multi)
	{
		Complex temp(0, multi);
		return temp;
	}
	double Re() const
	{
		return _re;
	}
	double Im() const
	{
		return _im;
	}
	double Re(const double real)
	{
		_re = real;
		return _re;
	}
	double Im(const double imaginary)
	{
		_im = imaginary;
		return _im;
	}
	double Modulus() const
	{
		return (sqrt(_re * _re + _im * _im));
	}
	double SquareModulus() const
	{
		return (_re * _re + _im * _im);
	}
	double Arg() const
	{
		return (atan2(_re, _im));
	}
	Complex Power(const double power) const
	{
		return (Power(*this, power));
	}
	Complex Root(const double root) const
	{
		return (Root(*this, root));
	}
	Complex Conjugate() const
	{
		Complex result = *this;
		result._im = -result._im;
		return result;
	}
	Complex& operator=(const Complex& complex)
	{
		this->_re = complex._re;
		this->_im = complex._im;
	}
	Complex operator+(const Complex& complex) const
	{
		Complex temp;
		temp._re = this->_re + complex._re;
		temp._im = this->_im + complex._im;
		return temp;
	}
	Complex operator-(const Complex& complex) const
	{
		Complex temp;
		temp._re = this->_re - complex._re;
		temp._im = this->_im - complex._im;
		return temp;
	}
	Complex operator*(const Complex& complex) const
	{
		Complex result;
		result._re = this->_re * complex._re - this->_im * complex._im;
		result._im = this->_re * complex._im + this->_im * complex._re;
		return result;
	}
	Complex operator/(const Complex& complex) const
	{
		Complex copy = complex;
		Complex result = *this * copy.Conjugate();
		double temp = copy.SquareModulus();
		result._re = result._re / temp;
		result._im = result._im / temp;
		return result;
	}
	Complex& operator=(const double real)
	{
		this->_re = real;
		this->_im = 0;
	}
	Complex operator+(const double real) const
	{
		Complex temp;
		temp._re = this->_re + real;
		return temp;
	}
	Complex operator-(const double real) const
	{
		Complex temp;
		temp._re = this->_re - real;
		return temp;
	}
	Complex operator*(const double real) const
	{
		Complex result;
		result._re = this->_re * real;
		result._im = this->_im * real;
		return result;
	}
	Complex operator/(const double real) const
	{
		Complex result = *this;
		result._re = result._re / real;
		result._im = result._im / real;
		return result;
	}
	Complex& operator-()
	{
		this->_re = -this->_re;
		this->_im = -this->_im;
		return (*this);
	}
	bool operator==(const Complex& complex) const
	{
		if (this->_re == complex._re && this->_im == complex._im) return true;
		else return false;
	}
	bool operator!=(const Complex& complex) const
	{
		if (this->_re != complex._re || this->_im != complex._im) return true;
		else return false;
	}
	bool operator==(const double real) const
	{
		if (this->_im == 0 && this->_re == real) return true;
		else return false;
	}
	bool operator!=(const double real) const
	{
		if (this->_im != 0 || this->_re != real) return true;
		else return false;
	}

	static double Re(const Complex& complex)
	{
		return complex._re;
	}
	static double Im(const Complex& complex)
	{
		return complex._im;
	}
	static double Modulus(const Complex& complex)
	{
		return (sqrt(complex._re * complex._re + complex._im * complex._im));
	}
	static double SquareModulus(const Complex& complex)
	{
		return (complex._re * complex._re + complex._im * complex._im);
	}
	static double Arg(const Complex& complex)
	{
		return (atan2(complex._re, complex._im));
	}
	static Complex ExpZ(const Complex& complex)
	{
		Complex result(cos(complex._im), sin(complex._im));
		result = result * pow(e, complex._re);
		return result;
	}
	static Complex ExpIZ(const Complex& complex)
	{
		Complex i(0, 1);
		Complex result = ExpZ(i * complex);
		return result;
	}
	static Complex Polar(const double r, const double theta)
	{
		Complex result(r * cos(AngleReduct(theta)), r * sin(AngleReduct(theta)));
		return result;
	}
	static Complex Power(const Complex& complex, const double power)
	{
		double r = Modulus(complex), theta = Arg(complex);
		r = pow(r, power), theta = AngleReduct(theta * power);
		Complex result(r, theta, true);
		return result;
	}
	static Complex Root(const Complex& complex, const double root)
	{
		double r = Modulus(complex), theta = Arg(complex);
		r = pow(r, 1.0 / root), theta = theta / root;
		Complex result(r, theta, true);
		return result;
	}
	static Complex Ln(const Complex& complex)
	{
		Complex result(log(Modulus(complex)), Arg(complex));
		return result;
	}
	static Complex Sin(const Complex& complex)
	{
		Complex copy = complex;
		Complex temp(0, 2);
		return ((ExpIZ(copy) - ExpIZ(-copy)) / (temp));
	}
	static Complex Cos(const Complex& complex)
	{
		Complex copy = complex;
		return ((ExpIZ(copy) + ExpIZ(-copy)) / 2.0);
	}
	static Complex Tan(const Complex& complex)
	{
		Complex copy = complex;
		Complex temp1 = ExpIZ(copy), temp2 = ExpIZ(-copy), unit(0, 1);
		Complex temp3 = (temp2 - temp1) * unit, temp4 = temp1 + temp2;
		return (temp3 / temp4);
	}
	static Complex Sinh(const Complex& complex)
	{
		Complex copy = complex;
		return ((ExpZ(copy) - ExpZ(-copy)) / 2.0);
	}
	static Complex Cosh(const Complex& complex)
	{
		Complex copy = complex;
		return ((ExpZ(copy) + ExpZ(-copy)) / 2.0);
	}

protected:
	double _re;
	double _im;
	static double AngleReduct(const double angle)
	{
		if (angle > -PI && angle <= PI) return angle;
		if (angle == -PI) return PI;
		int temp = angle / (2 * PI);
		double result = angle - temp * (2 * PI);
		if (result > PI) return (result - (2 * PI));
		return result;
	}

};
