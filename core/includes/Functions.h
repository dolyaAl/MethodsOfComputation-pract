#pragma once
#include <vector>
#include <cmath>
namespace Functions
{
	using ld = long double;
	class Function
	{
	public:
		Function() = default;
		/**
		 * Calculates the function value
		 * @param x: point at which value is calculated
		 * @return returns function value at x
		 */
		ld evaluate(ld);
		/**
		 * Counts function roots of odd multiplicity in range with precision parameter
		 * @param left: left edge of the range
		 * @param right: right edge of the range
		 * @param N: precision parameter (higher N - higher precision, default 10^3)
		 * @return returns roots count.
		 * If N = 0, returns empty vector
		*/
		int countRoots(ld, ld, int = pow(10, 3));
	};
	/**
	 * Separates range into N parts
	 * @param left: left edge of the range
	 * @param right: right edge of the range
	 * @param N: number of new ranges
	 * @return returns vector of smaller ranges.
	 * If N = 0, returns empty vector
	*/
	std::vector<std::pair<ld, ld>> separateRange(ld, ld, int);
}
