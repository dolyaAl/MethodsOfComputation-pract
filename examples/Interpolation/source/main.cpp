#include <iostream>
#include <iomanip>
#include <cmath>
#include "Functions.h"

using namespace std;

vector<pair<Functions::ld, Functions::ld>> table;
vector<pair<Functions::ld, Functions::ld>> nearest_table;
Functions::IFunction func;

Functions::ld f(Functions::ld x)
{
	return 1 - exp(-2 * x);
}

vector<pair<Functions::ld, Functions::ld>> getNearest(Functions::ld x, int N)
{
	vector<pair<Functions::ld, Functions::ld>> res;
	int i = 0;
	while(i < table.size() - 1 && table[i].second < x)
	{
		++i;
	}
	int k = i;
	int j = i + 1;
	while (res.size() < N && (j < table.size() && k >= 0))
	{
		if(abs(x - table[j].second) < abs(x - table[k].second))
		{
			res.push_back(table[j]);
			++j;
		}
		else
		{
			res.push_back(table[k]);
			--k;
		}
	}
	while(res.size() < N && j < table.size())
	{
		res.push_back(table[j]);
		++j;
	}
	while(res.size() < N && k >= 0)
	{
		res.push_back(table[k]);
		--k;
	}
	return res;
}

void doNewt(Functions::ld x, int N)
{
	func.initNewtonCoef(nearest_table, N);
	cout << "Newton: " << endl << "P(x) = " << func.evaluateNewtonInter(x) << endl;
	cout << " f(x) = " << f(x) << endl;
	cout << "|f(x) - P(x)| = " << abs(f(x) - func.evaluateNewtonInter(x)) << endl;
}

void doLagr(Functions::ld x, int N)
{
	func.initLagrangeCoef(nearest_table, N);
	cout << "Lagrange: " << endl << "P(x) = " << func.evaluateLagrangeInter(x) << endl;
	cout << " f(x) = " << f(x) << endl;
	cout << "|f(x) - P(x)| = " << abs(f(x) - func.evaluateLagrangeInter(x)) << endl;
}

void printNearestTable();
void doInterpolation()
{
	cout << "Enter x" << endl << ">>";
	Functions::ld x = 0;
	cin >> x;
	int N = 0;
	cout << "Enter degree " << "|||WARNING: DEGREE MUST BE <=" << table.size() - 1 << "|||" << endl << ">>";
	while (cin >> N)
	{
		if(N > table.size() - 1)
		{
			cout << "Incorrect polynomial degree, try again" << endl << ">>";
		}
		else
		{
			break;
		}
	}
	nearest_table = getNearest(x, N + 1);
	printNearestTable();
	system("pause");
	doNewt(x, N);
	doLagr(x, N);
}

void initTable()
{
	if(!table.empty())
	{
		table.clear();
	}
	Functions::ld A = 0;
	Functions::ld B = 0;
	Functions::ld x = 0;
	int m = 0;
	cout << "Enter M (units count)" << endl << ">>";
	cin >> m;
	cout << "Enter left edge" << endl << ">>";
	cin >> A;
	cout << "Enter right edge" << endl << ">>";
	cin >> B;
	--m;
	for (int i = 0; i <= m; ++i)
	{
		x = A + ((B - A) / m)*i;
		table.push_back({ f(x), x });
	}
}

void printTable()
{
	for (int i = 0; i < table.size(); ++i)
	{
		cout << "f(x_" << i << ") = " << table[i].first << " | x_" << i << " = " << table[i].second << endl;
	}
}

void printNearestTable()
{
	for (int i = 0; i < nearest_table.size(); ++i)
	{
		cout << "f(x_" << i << ") = " << nearest_table[i].first << " | x_" << i << " = " << nearest_table[i].second << endl;
	}
}

bool exit(bool exit = false)
{
	static bool ex;
	if(exit)
	{
		ex = true;
	}
	return ex;
}

void printMenu()
{
	cout << "The problem of algebraic interpolation. Var-5. f(x) = 1 - exp(-2x)" << endl;
	cout << "1. Interpolation"<<endl;
	cout << "2. Init table"<<endl;
	cout << "3. Print table"<<endl;
	cout << "4. Exit"<<endl;
}

void readCommand()
{
	int N = 0;
	cin >> N;
	switch(N)
	{
		case 1:
		{
			doInterpolation();
			system("pause");
			break;
		}
		case 2:
		{
			initTable();
			break;
		}
		case 3:
		{
			printTable();
			system("pause");
			break;
		}
		case 4:
		{
			exit(true);
			break;
		}
		default:
		{
			cout << "Incorrect command, try again";
		}
	}
}


int main()
{
	cout << setprecision(15);
	while(!exit())
	{
		printMenu();
		readCommand();
		system("cls");
	}
	return 0;
}
