#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "Functions.h"
#include <Matrix.hpp>

using namespace Functions;

Functions::Function::Function()
{
	fp = new FunctionParser_ld();
	ep = new Ev3::ExpressionParser();
	ep->SetVariableID("x", 0);
	fp->AddConstant("PI", std::acos(-1));
	fp->setEpsilon(1e-15);
}

Functions::Function::Function(const std::string& func):Function()
{
	int nerr = 0;
	derivates.push_back({ new Ev3::Expression(ep->Parse(func.c_str(), nerr)), fp });
	fp->Parse(func, "x");
}

Functions::Function::~Function()
{
	for (int i = 0; i < derivates.size(); ++i)
	{
		delete derivates[i].first;
		delete derivates[i].second;
	}
}

ld Function::f(ld x)
{
	return fp->Eval(&x);
}

ld Function::Df(ld x, int n)
{
	if(n < 0)
	{
		return INFINITY;
	}
	FunctionParser_ld* fup;
	Ev3::Expression* der;
	while (n >= derivates.size()) 
	{
		der = new Ev3::Expression(Ev3::Diff(*derivates.back().first, 0));
		fup = new FunctionParser_ld();	
		fup->AddConstant("PI", std::acos(-1));
		fup->Parse((*der)->ToString(), "x");
		fup->setEpsilon(1e-15);
		derivates.push_back({der, fup});
	}
	return derivates[n].second->Eval(&x);
}

void Functions::Function::updateFunc(const std::string& new_fun)
{
	int nerr = 0;
	for (int i = 1; i < derivates.size(); ++i)
	{
		delete derivates[i].first;
		delete derivates[i].second;
	}
	derivates.resize(1);
	delete derivates[0].first;
	derivates[0].first = new Ev3::Expression(ep->Parse(new_fun.c_str(), nerr));
	fp->Parse(new_fun, "x");
}

ld Functions::Function::findStartX(ld left, ld right)
{
	ld x = (left + right) / 2;
	for (int k = 2; k < 10000 && (f(x) * Df(x,2)) <= 0; ++k)
	{
		for (int j = 1; j < k; ++j)
		{
			x = left + (right - left) / k * j;
			if((f(x) * Df(x,2)) > 0)
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
		count += ((f(*i) * f(*(i + 1))) < 0.0);
		nulls += (f((*i)) == 0.0);
	}
	nulls += (f((*(ranges.end() - 1))) == 0);
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
		if(f(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		left = ranges[i];
		right = ranges[i + 1];
		if (f(left)*f(right) < 0)
		{
			count = 0;
			while(abs(left - right) >= 2*eps)
			{
				++count;
				mid = (left + right) / 2;
				if (f(left)*f(mid) < 0)
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
	if (f(ranges[ranges.size() - 1]) == 0)
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
		if(f(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		if (f(ranges[i+1])*f(ranges[i]) < 0)
		{
			cur = findStartX(ranges[i], ranges[i + 1]);
			prev = ranges[i];
			count = 0;
			while(abs(cur - prev) >= eps)
			{
				++count;
				prev = cur;
				cur = prev - f(prev) / Df(prev, 1);
			}
			roots.push_back(cur);
			logout << count << std::endl;
		}
	}
	if (f(ranges[ranges.size() - 1]) == 0)
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
		if(f(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		if (f(ranges[i+1])*f(ranges[i]) < 0)
		{
			cur = (ranges[i+1]+ranges[i])/2;
			x0 = findStartX(ranges[i], ranges[i+1]);
			prev = ranges[i];
			count = 0;
			while(abs(cur - prev) >= eps)
			{
				++count;
				prev = cur;
				cur = prev - f(prev) / Df(x0, 1);
			}
			roots.push_back(cur);
			logout << count << std::endl;
		}
	}
	if (f(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

std::vector<ld> Functions::Function::secantMethod(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	ld cur = 0;
	ld prev = 0;
	ld next = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(f(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		cur = (ranges[i+1]+ranges[i])/2;
		prev = ranges[i];
		if (f(ranges[i+1])*f(ranges[i]) < 0)
		{
			while(abs(cur - prev) >= eps)
			{
				next = cur - f(cur) / (f(cur) - f(prev)) * (cur - prev);
				prev = cur;
				cur = next;
			}
			roots.push_back(cur);
		}
	}
	if (f(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

ld Functions::Function::quadrfleftRect(ld left, ld right)
{
	return (right - left)*f(left);
}
ld Functions::Function::quadrfrightRect(ld left, ld right)
{
	return (right - left)*f(right);
}
ld Functions::Function::quadrfmidRect(ld left, ld right)
{
	return (right - left)*f((left + right)/2);
}
ld Functions::Function::quadrfTrapeze(ld left, ld right)
{
	return (right - left)*(f(left) + f(right))/2;
}
ld Functions::Function::quadrfSimpson(ld left, ld right)
{
	return (right - left)*(f(left) + 4*f((left+right)/2) + f(right)) / 6;
}
ld Functions::Function::quadrfThreeEighths(ld left, ld right)
{
	ld h = (right - left) / 3;
	return (right - left)*(f(left) + 3*f(left + h) + 3*f(left + 2*h) + f(right)) / 8;
}

ld Functions::Function::NIntleftRect(ld left, ld right, int M)
{
	ld h = (right - left) / M;
	ld res = 0;
	for (int i = 1; i <= M; ++i)
	{
		res += quadrfleftRect(left + (i - 1)*h, left + i * h);
	}
	return res;
}
ld Functions::Function::NIntrightRect(ld left, ld right, int M)
{
	ld h = (right - left) / M;
	ld res = 0;
	for (int i = 1; i <= M; ++i)
	{
		res += quadrfrightRect(left + (i - 1)*h, left + i*h);
	}
	return res;
}
ld Functions::Function::NIntmidRect(ld left, ld right, int M)
{
	ld h = (right - left) / M;
	ld res = 0;
	for (int i = 1; i <= M; ++i)
	{
		res += quadrfmidRect(left + (i - 1)*h, left + i*h);
	}
	return res;
}
ld Functions::Function::NIntTrapeze(ld left, ld right, int M)
{
	ld h = (right - left) / M;
	ld res = 0;
	for (int i = 1; i <= M; ++i)
	{
		res += quadrfTrapeze(left + (i - 1)*h, left + i*h);
	}
	return res;
}
ld Functions::Function::NIntSimpson(ld left, ld right, int M)
{
	ld h = (right - left) / M;
	ld res = 0;
	for (int i = 1; i <= M; ++i)
	{
		res += quadrfSimpson(left + (i - 1)*h, left + i*h);
	}
	return res;
}

std::vector<std::pair<ld, ld>> Functions::Function::initInterQF(std::vector<ld> weight_mom, std::vector<ld>& units)
{
	int N = units.size();
	Matrix<ld> matr = matrix.zeros<ld>(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			matr(i, j) = std::pow(units[j], i);
		}
	}
	std::vector<std::pair<ld, ld>> res;
	auto tmp = solveLin(matr, weight_mom, N);
	for (int i = 0; i < N; ++i)
	{
		res.push_back({ tmp[i], units[i]});
	}
	return res;
}

std::vector<ld> Functions::Function::findOrtPolynomCoefs(std::vector<ld> weight_mom)
{
	int N = weight_mom.size()/2;
	Matrix<ld> matr = matrix.zeros<ld>(N, N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			matr(i, j) = weight_mom[j+i];
		}
	}
	std::vector<ld> b(N);
	for (int j = N; j < 2*N; ++j)
	{
		b[j%N]= -weight_mom[j];
	}
	auto tmp = solveLin(matr, b, N);
	return tmp;
}

std::vector<ld> Functions::Function::OrtPolynomRoots(ld A, ld B, std::vector<ld> coefs)
{
	std::string polynom;
	std::stringstream str;
	str << std::setprecision(15);
	for(int i = 0; i < coefs.size(); ++i)
	{
		str << coefs[i];
		polynom += "(" + str.str() + ")*(x^" + std::to_string(i) + ")+";
		str.str(std::string());
	}
	polynom += "(x^" + std::to_string(coefs.size()) + ")";
	Function poly_f(polynom);
	return poly_f.bisectionMethod(1e-15, separateRange(A, B, (B - A) / (1e-3)));
}

std::vector<ld> Functions::Function::getGaussRoots(int n)
{
	static int N = -1;
	static std::vector<ld> roots = {};
	if(n == N)
	{
		return roots;
	}
	else
	{
		N = n;
	}
	Function legrPol(getNthLeg(n));
	roots = legrPol.secantMethod(1e-15, separateRange(-1, 1, 10000));
	return roots;
;
}

std::vector<ld> Functions::Function::getGaussCf(int n)
{
	static int N = -1;
	static std::vector<ld> cf;
	if(n == N)
	{
		return cf;
	}
	else
	{
		N = n;
	}
	auto lroots = getGaussRoots(n);
	Function tmpPol(getNthLeg(n-1));
	cf.clear();
	for (int i = 0; i < n; ++i)
	{
		cf.push_back(2.0*(1 - std::pow(lroots[i],2))/ std::pow(n*tmpPol.f(lroots[i]), 2));
	}
	return cf;
}

void Functions::Function::printGaussCf(std::vector<ld>& c, std::vector<ld>& t, ld A, ld B)
{
	ld cf = (B - A) / 2;
	for (int i = 0; i < c.size(); ++i)
	{
		std::cout << "C_" << std::to_string(i) << " = " << c[i] << 
			"\t A_" << std::to_string(i) << " = " << cf*c[i]<< std::endl;
		std::cout << "t_" << std::to_string(i) << " = " << t[i]<< 
			"\t x_" << std::to_string(i) << " = " << (cf*t[i]+(B+A)/2) << std::endl;
	}
}

ld Functions::Function::NIntInterQF(ld A, ld B, Weight& weight, int n, const std::vector<ld>& units)
{
	std::vector<ld> UN = units.empty() ? separateRange(A, B, n) : units;
	auto QF = initInterQF(weight.get_weight_moms(A, B,  n), UN);
	ld res = 0;
	for (int i = 0; i < QF.size(); ++i)
	{
		res += QF[i].first * f(QF[i].second);
	}
	return res;
}

ld Functions::Function::NIntMAGQF(ld A, ld B, Weight& weight, int n)
{
	auto units = OrtPolynomRoots(A, B, findOrtPolynomCoefs(weight.get_weight_moms(A, B, 2 * n)));
	return NIntInterQF(A, B, weight, n, units);
}

ld Functions::Function::NIntGauss(ld A, ld B, int n, bool print)
{
	auto lroots = getGaussRoots(n);
	std::vector<ld> C = getGaussCf(n);
	ld sum = 0;
	ld cf = (B - A) / 2;
	if (print)
	{
		printGaussCf(C, lroots, A, B);
	}
	for (int i = 0; i < n; ++i)
	{
		sum += C[i] * f(cf * lroots[i] + (B+A)/2);
	}
	return cf*sum;
}
ld Functions::Function::NIntMehler(int n)
{
	ld sum = 0;
	std::cout << "Meler cf = " << acos(-1) / n << std::endl;
	for (int i = 1; i <= n; ++i)
	{
		std::cout << "x_" << std::to_string(i) << " = " << std::cos((2 * i - 1) * std::acos(-1) / (2*n)) << std::endl;
		sum += f(std::cos((2 * i - 1) * std::acos(-1) / (2*n)));
	}
	return std::acos(-1)*sum/n;
}


ld Functions::IFunction::f(ld x)
{
	return evaluateNewtonInter(x);
}

void Functions::IFunction::initNewtonCoef(std::vector<std::pair<ld, ld>> table, int n)
{
	if(!m_NewtonInterpolationCoef.empty())
	{
		m_NewtonInterpolationCoef.clear();
		m_NewtonInterpolationUnits.clear();
	}
	for (int i = 0; i < n+1 && i < table.size(); ++i)
	{
		m_NewtonInterpolationUnits.push_back(table[i].second);
		m_NewtonInterpolationCoef.push_back(table[i].first - evaluateNewtonInter(table[i].second));
		for (int j = 0; j < i ; ++j)
		{
			m_NewtonInterpolationCoef[i] /= table[i].second - table[j].second;
		}
	}
}

ld Functions::IFunction::evaluateNewtonInter(ld x)
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

void Functions::IFunction::initLagrangeCoef(std::vector<std::pair<ld, ld>> table, int n)
{
	if(!m_LagrangeInterpolationCoef.empty())
	{
		m_LagrangeInterpolationCoef.clear();
		m_LagrangeInterpolationUnits.clear();
	}
	for (int i = 0; i < n+1 && i < table.size(); ++i)
	{
		m_LagrangeInterpolationUnits.push_back(table[i].second);
		m_LagrangeInterpolationCoef.push_back(table[i].first);
		for (int j = 0; j < n+1 && j < table.size(); ++j)
		{
			if(j != i)
			{
				m_LagrangeInterpolationCoef[i] /= table[i].second - table[j].second;
			}
		}
	}
}

ld Functions::IFunction::evaluateLagrangeInter(ld x)
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

std::vector<ld> Functions::IFunction::roots(ld eps, std::vector<ld>& ranges)
{
	std::vector<ld> roots;
	ld cur = 0;
	ld prev = 0;
	ld next = 0;
	for (int i = 0; i+1 < ranges.size(); ++i)
	{
		if(f(ranges[i]) == 0)
		{
			roots.push_back(ranges[i]);
			continue;
		}
		cur = (ranges[i+1]+ranges[i])/2;
		prev = ranges[i];
		if (f(ranges[i+1])*f(ranges[i]) < 0)
		{
			while(abs(cur - prev) >= eps)
			{
				next = cur - f(cur) / (f(cur) - f(prev)) * (cur - prev);
				prev = cur;
				cur = next;
			}
			roots.push_back(cur);
		}
	}
	if (f(ranges[ranges.size() - 1]) == 0)
	{
		roots.push_back(ranges[ranges.size() - 1]);
	}
	return roots;
}

std::vector<ld> Functions::Weight::get_weight_moms(ld A, ld B, int n)
{
	Function degr;
	for (int i = weight_moms.size(); i < n; ++i)
	{
		degr.updateFunc((*derivates[0].first)->ToString() + "*(" + "x^" + std::to_string(i) + ")");
		weight_moms.push_back(degr.NIntSimpson(A, B, 10000));
	}
	return std::vector<ld>(weight_moms.begin(), weight_moms.begin() + n);
}

std::vector<ld> Functions::separateRange(ld left, ld right, int N)
{
	std::vector<ld> res;
	N--;
	for (int i = 0; i <= N; ++i)
	{
		res.push_back(left + i * (right - left) / N);
	}
	return res;
}
long long int Functions::factorial(int n)
{
	long long int res = 1;
	for (int i = 1; i <= n; ++i)
	{
		res *= i;
	}
	return res;
}

std::vector<ld> Functions::solveLin(Matrix<ld> mat, std::vector<ld> b, int n)
{
	std::vector<ld> res(n);
	ld d = 0;
	for (int k = 0; k < n; k++)
	{
		for (int j = k + 1; j < n; j++)
		{
			d = mat.data_mat[j][k] / mat.data_mat[k][k];
			for (int i = k; i < n; i++)
			{
				mat.data_mat[j][i] = mat.data_mat[j][i] - d * mat.data_mat[k][i];
			}
			b[j] = b[j] - d * b[k];
		}
	}
	for (int k = n-1; k >= 0; k--)
	{
		d = 0;
		for (int j = k + 1; j < n; j++)
		{
			ld s = mat.data_mat[k][j] * res[j]; 
			d = d + s;
		}
		res[k] = (b[k] - d) / mat.data_mat[k][k];
	}
	return res;
}

std::string Functions::getNthLeg(int n)
{
	static std::vector<std::string> leg = { "1", "x" };
	for (int i = leg.size(); i <= n; ++i)
	{
		leg.push_back("(((2*" + std::to_string(i) + " - 1)/" + std::to_string(i) + ")*(" + leg[i - 1] + ")*x-("
			+ "((" + std::to_string(i) + " - 1)/" + std::to_string(i) + ")*(" + leg[i - 2] + ")))");
	}
	return leg[n];
}
