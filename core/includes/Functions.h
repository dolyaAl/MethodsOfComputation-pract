#pragma once
#include <vector>
#include <string>
#include <stack>
#include <cmath>
#include <expression.h>
#include <fparser.hh>
#include <parser.h>
#include <Matrix.hpp>
namespace Functions
{
	using ld = long double;
	class Weight;
	class Function
	{
	protected:
		FunctionParser_ld* fp;
		Ev3::ExpressionParser* ep;
		std::vector< std::pair<Ev3::Expression*, FunctionParser_ld* >> derivates;
		
	public:
		Function();
		Function(const std::string& fu);
		~Function();
		/**
		 * Calculates the function value
		 * @param x: point at which value is calculated
		 * @return returns function value at x
		 */
		ld f(ld x);
		/**
		 * Calculates the function nth derivative value
		 * @param x: point at which value is calculated
		 * @return returns derivative value at x
		 */
		ld Df(ld x, int n);
		void updateFunc(const std::string& new_fun);
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
		 * @param use_iner_polynom: set true to use Newton's interpolation polynom for evaluating function
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> secantMethod(ld, std::vector<ld>&);

		ld quadrfleftRect(ld left, ld right);
		ld quadrfrightRect(ld left, ld right);
		ld quadrfmidRect(ld left, ld right);
		ld quadrfTrapeze(ld left, ld right);
		ld quadrfSimpson(ld left, ld right);
		ld quadrfThreeEighths(ld left, ld right);

		ld NIntleftRect(ld left, ld right, int M);
		ld NIntrightRect(ld left, ld right, int M);
		ld NIntmidRect(ld left, ld right, int M);
		ld NIntTrapeze(ld left, ld right, int M);
		ld NIntSimpson(ld left, ld right, int M);
		ld RungeAccSimpson(ld j1, ld j2, int L);
	private: 
		std::vector<std::pair<ld, ld>> initInterQF(std::vector<ld> weight_mom, std::vector<ld>& units);
		std::vector<ld> findOrtPolynomCoefs(std::vector<ld> weight_mom);
		std::vector<ld> OrtPolynomRoots(ld A, ld B, std::vector<ld> coefs);
		std::vector<ld> getGaussRoots(int n);
		std::vector<ld> getGaussCf(int n);
		void printGaussCf(std::vector<ld>& c, std::vector<ld>& t, ld A, ld B);
	public:
		ld NIntInterQF(ld A, ld B, Weight& weight, int n, const std::vector<ld>& units = std::vector<ld>());
		ld NIntMAGQF(ld A, ld B, Weight& weight, int n);
		ld NIntGauss(ld A, ld B, int n, bool print = false);
		ld NIntMeler(int n);


		
	};
	class IFunction
	{
		std::vector<ld> m_NewtonInterpolationCoef;
		std::vector<ld> m_NewtonInterpolationUnits;
		std::vector<ld> m_LagrangeInterpolationCoef;
		std::vector<ld> m_LagrangeInterpolationUnits;
	public:
		/**
		 * Calculates the function value
		 * @param x: point at which value is calculated
		 * @return returns function value at x
		 */
		ld f(ld x);

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
		/**
		 * Finds function roots of odd multiplicity
		 * @param eps: sets the accuracy of the calculation
		 * @param ranges: vector of the partition points of the interval on which to look for roots.
		 * @param use_iner_polynom: set true to use Newton's interpolation polynom for evaluating function
		 * get it from separateRange function in Functions namespace
		 * @return vector of roots
		*/
		std::vector<ld> roots(ld, std::vector<ld>&);
	};
	class Weight : public Function
	{
		std::vector<ld> weight_moms;
	public:
		Weight(const std::string weight) :Function(weight) {};
		std::vector<ld> get_weight_moms(ld A, ld B, int n);
	};
	/*
	* 
	 * Separates range into N parts
	 * @param left: left edge of the range
	 * @param right: right edge of the range
	 * @param N: number of new ranges
	 * @return returns vector of smaller ranges.
	 * If N = 0, returns empty vector
	*/
	std::vector<ld> separateRange(ld left, ld right, int N);
	long long int factorial(int n);
	std::vector<ld> solveLin(Matrix<ld> mat, std::vector<ld> b, int n);
	std::string getNthLeg(int n);
}
