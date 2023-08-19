#pragma once
#include <future>
#include <string>
#include <thread>
#include <vector>
#include "BIGINTEGER.h"
#include "COMPLEX.h"
#include "LINEAR.h"

class Fast
{
public:
	static BigInteger A(const long long n, const long long r)
	{
		if (n < 0 || r < 0) return 0;
		if (n >= r) return (Factorial(n - r + 1, n));
		if (n <= r) return (Factorial(r - n + 1, r));
	}
	static BigInteger C(const long long n, const long long r)
	{
		if (n < 0 || r < 0) return 0;
		if (n >= r) {
			if (r > n / 2) return (C(n, n - r));
			else return (Factorial(n - r + 1, n) / Factorial(r));
		}
		if (n <= r) {
			if (n > r / 2) return (C(r, r - n));
			else return (Factorial(r - n + 1, r) / Factorial(n));
		}
	}
	static BigInteger Factorial(const long long num)
	{
		BigInteger result;
		if (num < 500) result = SimFac(1, num, 1);
		else {
			long long unit = num / 8;
			std::future<BigInteger> f1 = std::async(std::launch::async, SimFac, 1, unit, 1);
			std::future<BigInteger> f2 = std::async(std::launch::async, SimFac, unit + 1, 2 * unit, 1);
			std::future<BigInteger> f3 = std::async(std::launch::async, SimFac, 2 * unit + 1, 3 * unit, 1);
			std::future<BigInteger> f4 = std::async(std::launch::async, SimFac, 3 * unit + 1, 4 * unit, 1);
			std::future<BigInteger> f5 = std::async(std::launch::async, SimFac, 4 * unit + 1, 5 * unit, 1);
			std::future<BigInteger> f6 = std::async(std::launch::async, SimFac, 5 * unit + 1, 6 * unit, 1);
			std::future<BigInteger> f7 = std::async(std::launch::async, SimFac, 6 * unit + 1, 7 * unit, 1);
			std::future<BigInteger> f8 = std::async(std::launch::async, SimFac, 7 * unit + 1, num, 1);
			result = f1.get() * f2.get() * f3.get() * f4.get() * f5.get() * f6.get() * f7.get() * f8.get();
		}
		return result;
	}
	static BigInteger DoubleFactorial(const long long num)
	{
		BigInteger result;
		if (num < 1000) {
			if (num % 2 == 0) result = SimFac(2, num, 2); //若是偶数
			else result = SimFac(1, num, 2); //若是奇数
		}
		else {
			if (num % 2 == 0) { //若是偶数
				long long unit = num / 8;
				if (unit % 2 != 0) unit = unit + 1; //保证unit是偶数
				std::future<BigInteger> f1 = std::async(std::launch::async, SimFac, 2, unit, 2);
				std::future<BigInteger> f2 = std::async(std::launch::async, SimFac, unit + 2, 2 * unit, 2);
				std::future<BigInteger> f3 = std::async(std::launch::async, SimFac, 2 * unit + 2, 3 * unit, 2);
				std::future<BigInteger> f4 = std::async(std::launch::async, SimFac, 3 * unit + 2, 4 * unit, 2);
				std::future<BigInteger> f5 = std::async(std::launch::async, SimFac, 4 * unit + 2, 5 * unit, 2);
				std::future<BigInteger> f6 = std::async(std::launch::async, SimFac, 5 * unit + 2, 6 * unit, 2);
				std::future<BigInteger> f7 = std::async(std::launch::async, SimFac, 6 * unit + 2, 7 * unit, 2);
				std::future<BigInteger> f8 = std::async(std::launch::async, SimFac, 7 * unit + 2, num, 2);
				result = f1.get() * f2.get() * f3.get() * f4.get() * f5.get() * f6.get() * f7.get() * f8.get();
			}
			else { //若是奇数
				long long unit = num / 8;
				if (unit % 2 == 0) unit = unit + 1; //保证unit是奇数
				std::future<BigInteger> f1 = std::async(std::launch::async, SimFac, 1, unit, 2);
				std::future<BigInteger> f2 = std::async(std::launch::async, SimFac, unit + 2, 2 * unit - 1, 2);
				std::future<BigInteger> f3 = std::async(std::launch::async, SimFac, 2 * unit + 1, 3 * unit, 2);
				std::future<BigInteger> f4 = std::async(std::launch::async, SimFac, 3 * unit + 2, 4 * unit - 1, 2);
				std::future<BigInteger> f5 = std::async(std::launch::async, SimFac, 4 * unit + 1, 5 * unit, 2);
				std::future<BigInteger> f6 = std::async(std::launch::async, SimFac, 5 * unit + 2, 6 * unit - 1, 2);
				std::future<BigInteger> f7 = std::async(std::launch::async, SimFac, 6 * unit + 1, 7 * unit, 2);
				std::future<BigInteger> f8 = std::async(std::launch::async, SimFac, 7 * unit + 2, num, 2);
				result = f1.get() * f2.get() * f3.get() * f4.get() * f5.get() * f6.get() * f7.get() * f8.get();
			}
		}
		return result;
	}
	static BigInteger Fibonacci(const long long num)
	{
		Linear::Matrix<BigInteger> coeff(2);
		coeff.SetElmn(1, 1, 1), coeff.SetElmn(1, 1, 2), coeff.SetElmn(1, 2, 1), coeff.SetElmn(0, 2, 2);
		coeff = MatrixPower(coeff, num);
		return (coeff.GetElmn(2, 2));
	}
	static BigInteger ToBase10(const std::string& num, const int base)
	{
		int len = num.length();
		BigInteger bbase(base), result(0);
		for (int i = len - 1; i >= 0; i--) {
			if ((char)num[i] >= '0' && (char)num[i] <= '9') result = result + bbase.Power(len - i - 1) * (long long)((char)num[i] - '0');
			if ((char)num[i] >= 'A' && (char)num[i] <= 'Z') result = result + bbase.Power(len - i - 1) * (long long)((char)num[i] - 'A' + 10);
		}
		return result;
	}
	static std::string ToBase(const BigInteger& num, const int base)
	{
		char unit[37] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		std::string result;
		BigInteger copy(num), bbase(base), remain(0);
		while (copy != (long long)0) {
			copy = BigInteger::Division(copy, bbase, remain);
			result.insert(result.begin(), unit[(long long)remain]);
		}
		return result;
	}

	template <typename T>
	static T Power(const T& back, const long long power)
	{
		T result = 1, copy = back;
		if (power > 0) {
			long long pow = power;
			while (pow) {
				if (pow & 1) result = result * copy;
				copy = copy * copy;
				pow = pow / 2;
			}
			return result;
		}
		else return back;
	}
	template <typename T>
	static Linear::Matrix<T> MatrixPower(const Linear::Matrix<T>& back, const long long power)
	{
		Linear::Matrix<T> result = back, copy = back;
		result.UnitMatrix(back.Rank());
		if (power > 0) {
			long long pow = power;
			while (pow) {
				if (pow & 1) result = result * copy;
				copy = copy * copy;
				pow = pow / 2;
			}
			return result;
		}
		else return back;
	}


protected:
	static BigInteger SimFac(const long long start, const long long end, int step)
	{
		BigInteger result(1);
		for (long long i = start; i <= end; i = i + step) result = result * i;
		return result;
	}
	static BigInteger Factorial(const long long start, const long long end)
	{
		BigInteger result;
		if (end - start < 500) result = SimFac(start, end, 1);
		else {
			long long unit = end - start / 8;
			std::future<BigInteger> f1 = std::async(std::launch::async, SimFac, start, start + unit, 1);
			std::future<BigInteger> f2 = std::async(std::launch::async, SimFac, start + unit + 1, start + 2 * unit, 1);
			std::future<BigInteger> f3 = std::async(std::launch::async, SimFac, start + 2 * unit + 1, start + 3 * unit, 1);
			std::future<BigInteger> f4 = std::async(std::launch::async, SimFac, start + 3 * unit + 1, start + 4 * unit, 1);
			std::future<BigInteger> f5 = std::async(std::launch::async, SimFac, start + 4 * unit + 1, start + 5 * unit, 1);
			std::future<BigInteger> f6 = std::async(std::launch::async, SimFac, start + 5 * unit + 1, start + 6 * unit, 1);
			std::future<BigInteger> f7 = std::async(std::launch::async, SimFac, start + 6 * unit + 1, start + 7 * unit, 1);
			std::future<BigInteger> f8 = std::async(std::launch::async, SimFac, start + 7 * unit + 1, end, 1);
			result = f1.get() * f2.get() * f3.get() * f4.get() * f5.get() * f6.get() * f7.get() * f8.get();
		}
		return result;
	}
	static void Multiply(const BigInteger A[3][3], const BigInteger B[3][3], BigInteger C[3][3])
	{
		C[1][1] = A[1][1] * B[1][1] + A[1][2] * B[2][1];
		C[1][2] = A[1][1] * B[1][2] + A[1][2] * B[2][2];
		C[2][1] = A[2][1] * B[1][1] + A[2][2] * B[2][1];
		C[2][2] = A[2][1] * B[1][2] + A[2][2] * B[2][2];
	}

};
