#include "Functions.h"
using namespace Functions;
ld Function::evaluate(ld x)
{
	return x * sin(x) - 1;
}

int Functions::Function::countRoots(ld left, ld right, int N)
{
	std::vector <std::pair<ld, ld>> ranges = separateRange(left, right, N);
	int count = 0;
	int nulls = 0;
	for (auto i = ranges.begin(); i != ranges.end(); ++i)
	{
		count += (evaluate((*i).first) * evaluate((*i).second) < 0.0);
		nulls += (evaluate((*i).first) == 0.0);
	}
	nulls += (evaluate((*(ranges.end() - 1)).second) == 0);
	return count + nulls;
}

std::vector<std::pair<ld, ld>> Functions::separateRange(ld left, ld right, int N)
{
	std::vector<std::pair<ld, ld>> res;
	ld step = (right - left) / N;
	for (int i = 1; i < N; ++i)
	{
		res.push_back({ left + (i - 1) * step, left + i * step });
	}
	return res;
}
