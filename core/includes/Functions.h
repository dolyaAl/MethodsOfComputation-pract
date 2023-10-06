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
		/**
		* Constructs function by string fu
		* @param fu: symbolic representation of the function
		*/
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
		/**
		* Sets new symbolic representation to the function
		* @param new_fun: new symbolic representation of the function
		*/
		void updateFunc(const std::string& new_fun);
		/**
		 * Counts function roots of odd multiplicity in range with precision parameter
		 * @param left: left edge of the range
		 * @param right: right edge of the range
		 * @param N: precision parameter (higher N - higher precision, default 10^3)
		 * @return returns roots count.
		 * If N = 0, returns empty vector
		*/

//function roots

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

//quadratic formulas

		/**
		* Calculates the left rectangle formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfleftRect(ld left, ld right);
		/**
		* Calculates the rigth rectangle formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfrightRect(ld left, ld right);
		/**
		* Calculates the mid rectangle formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfmidRect(ld left, ld right);
		/**
		* Calculates the trapeze formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfTrapeze(ld left, ld right);
		/**
		* Calculates the Simpson's formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfSimpson(ld left, ld right);
		/**
		* Calculates the Simpson's three eights formula for the function on the interval [left, right]
		* @param left: left edge
		* @param right: right edge
		*/
		ld quadrfThreeEighths(ld left, ld right);

//Integration using quadratic formulas


		/**
		* Calculates the function integral on the interval [left, right] by left rectangle method
		* @param left: left edge
		* @param right: right edge
		*/
		ld NIntleftRect(ld left, ld right, int M);
		/**
		* Calculates the function integral on the interval [left, right] by right rectangle method
		* @param left: left edge
		* @param right: right edge
		*/
		ld NIntrightRect(ld left, ld right, int M);
		/**
		* Calculates the function integral on the interval [left, right] by mid rectangle method
		* @param left: left edge
		* @param right: right edge
		*/
		ld NIntmidRect(ld left, ld right, int M);
		/**
		* Calculates the function integral on the interval [left, right] by trapeze method
		* @param left: left edge
		* @param right: right edge
		*/
		ld NIntTrapeze(ld left, ld right, int M);
		/**
		* Calculates the function integral on the interval [left, right] by Simpson's method
		* @param left: left edge
		* @param right: right edge
		*/
		ld NIntSimpson(ld left, ld right, int M);

	private: 
		std::vector<std::pair<ld, ld>> initInterQF(std::vector<ld> weight_mom, std::vector<ld>& units);
		std::vector<ld> findOrtPolynomCoefs(std::vector<ld> weight_mom);
		std::vector<ld> OrtPolynomRoots(ld A, ld B, std::vector<ld> coefs);
		std::vector<ld> getGaussRoots(int n);
		std::vector<ld> getGaussCf(int n);
		void printGaussCf(std::vector<ld>& c, std::vector<ld>& t, ld A, ld B);
	public:
		/**
		* Calculates the function integral with weight on the interval [A, B] by interpolation quadratic formula 
		* @param A: left edge
		* @param B: right edge
		* @param weight: weight
		* @param n: the number units of the separation of the interval [A, B]
		* @param units: units of the separation of interval [A,B] (by deafault an equidistant separation into n units is used)
		*/
		ld NIntInterQF(ld A, ld B, Weight& weight, int n, const std::vector<ld>& units = std::vector<ld>());
		/**
		* Calculates the function integral with weight on the interval [A, B] by quadratic formula of maximum algebraic accuracy degree 
		* @param A: left edge
		* @param B: right edge
		* @param weight: weight
		* @param n: algebraic accuracy degree 2n-1
		*/
		ld NIntMAGQF(ld A, ld B, Weight& weight, int n);
		/**
		* Calculates the function integral with weight on the interval [A, B] by Gaussian quadratic formula 
		* @param A: left edge
		* @param B: right edge
		* @param n: algebraic accuracy degree 2n-1
		* @param print: if true prints Gaussian coefficients
		*/
		ld NIntGauss(ld A, ld B, int n, bool print = false);
		/**
		* Calculates the function integral with weight on the interval [A, B] by Mehler quadratic formula 
		* @param n: algebraic accuracy degree 2n-1
		*/
		ld NIntMehler(int n);


		
	};

//Interpolated function

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

//Weight

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
	/**
	* returns Nth degree Legendre polynomial symbolic representation
	* @param n: polynomial degree
	*/
	std::string getNthLeg(int n);
}
