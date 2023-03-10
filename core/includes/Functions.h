#pragma once
#include <vector>
#include <cmath>
namespace Functions
{
	using ld = long double;
	class Function
	{
		std::vector<ld> m_NewtonInterpolationCoef;
		std::vector<ld> m_NewtonInterpolationUnits;
		std::vector<ld> m_LagrangeInterpolationCoef;
		std::vector<ld> m_LagrangeInterpolationUnits;
	public:
		Function() = default;
		/**
		 * Calculates the function value
		 * @param x: point at which value is calculated
		 * @return returns function value at x
		 */
		ld evaluate(ld x);
		/**
		 * Calculates the function derivative value
		 * @param x: point at which value is calculated
		 * @return returns derivative value at x
		 */
		ld evaluateDer(ld x);
		/**
		 * Calculates the function second derivative value
		 * @param x: point at which value is calculated
		 * @return returns derivative value at x
		 */
		ld evaluateSecondDer(ld x);
		/**
		 * Counts function roots of odd multiplicity in range with precision parameter
		 * @param left: left edge of the range
		 * @param right: right edge of the range
		 * @param N: precision parameter (higher N - higher precision, default 10^3)
		 * @return returns roots count.
		 * If N = 0, returns empty vector
		*/
		int countRoots(ld left, ld right, int = pow(10, 3));
		/**
		 * Finds function roots of odd multiplicity using bisection method
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> bisectionMethod(ld eps, std::vector<ld>& ranges);
	private:
		/**
		 * Auxiliary function for newton's method
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		ld findStartX(ld left, ld right);
	public:
		/**
		 * Finds function roots of odd multiplicity using newton's method
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> newtonMethod(ld, std::vector<ld>&);
		/**
		 * Finds function roots of odd multiplicity using newton's modificated method
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> newtonModMethod(ld, std::vector<ld>&);

		/**
		 * Finds function roots of odd multiplicity using method
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> secantMethod(ld, std::vector<ld>&);
		/**
		 * Initializes coefficients for Newton's interpolation polynomial by table of values |f(x_i)|x_i|
		 * @param table: table[i].first - func value, table[i].second - point
		 * @param n: polynomial degree
		*/
		void initNewtonCoef(std::vector<std::pair<ld, ld>> table, int n);
		/**
		 * Evaluates function value by Newton's polynomial in point x
		 * note: initialize coefficients before using, use initNewtonCoef function
		 * @param x: point at which value is evaluated
		*/
		ld evaluateNewtonInter(ld x);
		/**
		 * Initializes coefficients for Lagrange interpolation polynomial by table of values |f(x_i)|x_i|
		 * @param table: table[i].first - func value, table[i].second - point
		 * @param n: polynomial degree
		*/
		void initLagrangeCoef(std::vector<std::pair<ld, ld>> table, int n);
		/**
		 * Evaluates function value by Lagrange polynomial in point x
		 * note: initialize coefficients before using, use initLagrangeCoef function
		 * @param x: point at which value is evaluated
		*/
		ld evaluateLagrangeInter(ld x);
	};
	/**
	 * Separates range into N parts
	 * @param left: left edge of the range
	 * @param right: right edge of the range
	 * @param N: number of new ranges
	 * @return returns vector of smaller ranges.
	 * If N = 0, returns empty vector
	*/
	std::vector<ld> separateRange(ld left, ld right, int N);
}
