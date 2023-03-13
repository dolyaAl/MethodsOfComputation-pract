#include <fstream>
#include <iostream>
#include "Functions.h"
using namespace Functions;
ld Function::evaluate(ld x)
{
	return x * sin(x) - 1;
}

ld Functions::Function::evaluateDer(ld x)
{
	constexpr ld h = 1e-10;
	return (evaluate(x+h) - evaluate(x - h))/(2*h);
}

ld Functions::Function::evaluateSecondDer(ld x)
{
	constexpr ld h = 1e-10;
	return (evaluateDer(x+h) - evaluateDer(x - h))/(2*h);
}

ld Functions::Function::findStartX(ld left, ld right)
{
	ld x = (left + right) / 2;
	for (int k = 2; k < 10000 && (evaluate(x) * evaluateSecondDer(x)) <= 0; ++k)
	{
		for (int j = 1; j < k; ++j)
		{
			x = left + (right - left) / k * j;
			if((evaluate(x) * evaluateSecondDer(x)) > 0)
			{
				break;
			}
		}
		if(k == 9999)
		{
			std::cout << "Something wrong check function findStartX" << std::endl;
		}
	}
	return x;
}

int Functions::Function::countRoots(ld left, ld right, int N)
{
	std::vector <ld> ranges = separateRange(left, right, N);
	int count = 0;
	int nulls = 0;
	for (auto i = ranges.begin(); i + 1 != ranges.end(); ++i)
	{
		count += ((evaluate(*i) * evaluate(*(i + 1))) < 0.0);
		nulls += (evaluate((*i)) == 0.0);
	}
	nulls += (evaluate((*(ranges.end() - 1))) == 0);
	return count + nulls;
}

std::vector<ld> Functions::Function::bisectionMethod(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	std::ofstream logout("bisection.txt");
	int count = 0;
	ld left = 0;
	ld right = 0;
	ld mid = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(evaluate(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		left = ranges[i];
		right = ranges[i + 1];
		if (evaluate(left)*evaluate(right) < 0)
		{
			count = 0;
			while(abs(left - right) >= 2*eps)
			{
				++count;
				mid = (left + right) / 2;
				if (evaluate(left)*evaluate(mid) < 0)
				{
					right = mid;
				}
				else
				{
					left = mid;
				}
			}
			roots.push_back((left + right) / 2);
			logout << count << std::endl;
		}
	}
	if (evaluate(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

std::vector<ld> Functions::Function::newtonMethod(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	std::ofstream logout("newton.txt");
	int count = 0;
	ld cur = 0;
	ld prev = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(evaluate(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		if (evaluate(ranges[i+1])*evaluate(ranges[i]) < 0)
		{
			cur = findStartX(ranges[i], ranges[i + 1]);
			prev = ranges[i];
			count = 0;
			while(abs(cur - prev) >= eps)
			{
				++count;
				prev = cur;
				cur = prev - evaluate(prev) / evaluateDer(prev);
			}
			roots.push_back(cur);
			logout << count << std::endl;
		}
	}
	if (evaluate(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

std::vector<ld> Functions::Function::newtonModMethod(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	std::ofstream logout("newtonMod.txt");
	int count = 0;
	ld cur = 0;
	ld prev = 0;
	ld x0 = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(evaluate(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		if (evaluate(ranges[i+1])*evaluate(ranges[i]) < 0)
		{
			cur = (ranges[i+1]+ranges[i])/2;
			x0 = findStartX(ranges[i], ranges[i+1]);
			prev = ranges[i];
			count = 0;
			while(abs(cur - prev) >= eps)
			{
				++count;
				prev = cur;
				cur = prev - evaluate(prev) / evaluateDer(x0);
			}
			roots.push_back(cur);
			logout << count << std::endl;
		}
	}
	if (evaluate(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

std::vector<ld> Functions::Function::secantMethod(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	std::ofstream logout("secant.txt");
	int count = 0;
	ld cur = 0;
	ld prev = 0;
	ld next = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(evaluate(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		cur = (ranges[i+1]+ranges[i])/2;
		prev = ranges[i];
		if (evaluate(ranges[i+1])*evaluate(ranges[i]) < 0)
		{
			count = 0;
			while(abs(cur - prev) >= eps)
			{
				++count;
				next = cur - evaluate(cur) / (evaluate(cur) - evaluate(prev)) * (cur - prev);
				prev = cur;
				cur = next;
			}
			roots.push_back(cur);
			logout << count << std::endl;
		}
	}
	if (evaluate(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

void Functions::Function::initNewtonCoef(std::vector<std::pair<ld, ld>> table, int n)
{
	if(!m_NewtonInterpolationCoef.empty())
	{
		m_NewtonInterpolationCoef.clear();
		m_NewtonInterpolationUnits.clear();
	}
	for (int i = 0; i < n && i < table.size(); ++i)
	{
		m_NewtonInterpolationUnits.push_back(table[i].second);
		m_NewtonInterpolationCoef.push_back(table[i].first - evaluateNewtonInter(table[i].second));
		for (int j = 0; j < i ; ++j)
		{
			m_NewtonInterpolationCoef[i] /= table[i].second - table[j].second;
		}
	}
}

ld Functions::Function::evaluateNewtonInter(ld x)
{
	ld res = 0;
	ld tmp = 0;
	for (int i = 0; i < m_NewtonInterpolationCoef.size(); ++i)
	{
		tmp = m_NewtonInterpolationCoef[i];
		for (int j = 0; j < i ; ++j)
		{
			tmp *= x - m_NewtonInterpolationUnits[j];
		}
		res += tmp;
	}
	return res;
}

void Functions::Function::initLagrangeCoef(std::vector<std::pair<ld, ld>> table, int n)
{
	if(!m_LagrangeInterpolationCoef.empty())
	{
		m_LagrangeInterpolationCoef.clear();
		m_LagrangeInterpolationUnits.clear();
	}
	for (int i = 0; i < n && i < table.size(); ++i)
	{
		m_LagrangeInterpolationUnits.push_back(table[i].second);
		m_LagrangeInterpolationCoef.push_back(table[i].first);
		for (int j = 0; j < n && j < table.size(); ++j)
		{
			if(j != i)
			{
				m_LagrangeInterpolationCoef[i] /= table[i].second - table[j].second;
			}
		}
	}
}

ld Functions::Function::evaluateLagrangeInter(ld x)
{
	ld res = 0;
	ld tmp = 0;
	for (int i = 0; i < m_LagrangeInterpolationUnits.size(); ++i)
	{
		tmp = m_LagrangeInterpolationCoef[i];
		for(int j = 0; j < m_LagrangeInterpolationUnits.size(); ++j)
		{
			if (j != i)
			{
				tmp *= x - m_LagrangeInterpolationUnits[j];
			}
		}
		res += tmp;
	}
	return res;
}

std::vector<ld> Functions::separateRange(ld left, ld right, int N)
{
	std::vector<ld> res;
	for (int i = 1; i <= N; ++i)
	{
		res.push_back(left + i * (right - left) / N);
	}
	return res;
}
