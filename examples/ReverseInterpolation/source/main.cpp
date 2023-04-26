#include <iostream>
#include <iomanip>
#include <cmath>
#include "Functions.h"

using namespace std;

vector<pair<Functions::ld, Functions::ld>> table;
Functions::IFunction func;
void printTable(vector<pair<Functions::ld, Functions::ld>>);

Functions::ld f(Functions::ld x)
{
	return sin(5*acos(-1)*x);
}

vector<pair<Functions::ld, Functions::ld>> getNearest(vector<pair<Functions::ld, Functions::ld>> itable, Functions::ld x, int N)
{
	vector<pair<Functions::ld, Functions::ld>> res;
	int i = 0;
	int min_delta_indx = 0;
	Functions::ld delta_min = INFINITY;
	Functions::ld delta = 0;
	while(i < itable.size() - 1)
	{
		delta = abs(itable[i].second - x);
		if(delta < delta_min)
		{
			delta_min = delta;
			min_delta_indx = i;
		}
		++i;
	}
	int k = min_delta_indx;
	int j = min_delta_indx + 1;
	while (res.size() < N+1 && (j < itable.size() && k >= 0))
	{
		if(abs(x - itable[j].second) < abs(x - itable[k].second))
		{
			res.push_back(itable[j]);
			++j;
		}
		else
		{
			res.push_back(itable[k]);
			--k;
		}
	}
	while(res.size() < N+1 && j < itable.size())
	{
		res.push_back(itable[j]);
		++j;
	}
	while(res.size() < N+1 && k >= 0)
	{
		res.push_back(itable[k]);
		--k;
	}
	return res;
}
vector<pair<Functions::ld, Functions::ld>> shiftTable(vector<pair<Functions::ld, Functions::ld>> tab, Functions::ld F)
{
	vector<pair<Functions::ld, Functions::ld>> res;
	for(int i = 0; i < tab.size(); ++i)
	{
		res.push_back({ tab[i].first - F, tab[i].second });
	}
	return res;
}
vector<pair<Functions::ld, Functions::ld>> reverseTable(vector<pair<Functions::ld, Functions::ld>> tab)
{
	vector<pair<Functions::ld, Functions::ld>> res;
	for(int i = 0; i < tab.size(); ++i)
	{
		res.push_back({ tab[i].second, tab[i].first});
	}
	return res;
}


vector<vector<pair<Functions::ld, Functions::ld>>> getMonotonTable(vector<pair<Functions::ld, Functions::ld>> tab)
{
	vector<vector<pair<Functions::ld, Functions::ld>>> res;
	vector<pair<Functions::ld, Functions::ld>> temp;
	temp.push_back(tab.front());
	for(int i = 1; i < tab.size(); ++i)
	{
		if (temp.back() < tab[i])
		{
			while(i < tab.size() && temp.back() < tab[i])
			{
				temp.push_back(tab[i]);
				++i;
			}
		}
		else
		{
			while(i < tab.size() && temp.back() >= tab[i])
			{
				temp.push_back(tab[i]);
				++i;
			}
		}
		if(!temp.empty())
		{
			res.push_back(temp);
			temp.clear();
			if(i < tab.size())
			{
				temp.push_back(tab[i]);
			}
		}
	}
	return res;
}

void doFirstMethod(Functions::ld F, int N)
{
	cout << "First method--------------------------------------------------------------------" << endl;
	auto tables = getMonotonTable(table);
	for (int i = 0; i < tables.size(); ++i)
	{
		if(!tables[i].empty())
		{
			auto reverse_t = reverseTable(tables[i]);
			//printTable(tables[i]);
			func.initNewtonCoef(getNearest(reverse_t, F, N), N);
			cout << "P(F) = " << func.evaluateNewtonInter(F) << endl;
			cout << "|f(P(F)) - F| = " << abs(f(func.evaluateNewtonInter(F)) - F) << endl;
		}
	}
	cout << endl;
	cout << endl;
}
/*
void doFirstMethod(Functions::ld F, int N)
{
	auto reverse_t = reverseTable(table);
	func.initNewtonCoef(getNearest(reverse_t, F, N), N);
	cout << "First method" << endl << "P(F) = " << func.evaluateNewtonInter(F) << endl;
	cout << "|f(P(F)) - F| = " << abs(f(func.evaluateNewtonInter(F)) - F) << endl;
}*/

void doSecondMethod(Functions::ld A, Functions::ld B,Functions::ld eps, int N, Functions::ld F)
{
	cout << "Second method--------------------------------------------------------------------" << endl;
	auto monoton_tabs = getMonotonTable(table);
	for (auto& tab : monoton_tabs)
	{
		if ((tab.front().first <= F && tab.back().first >= F) || (tab.front().first >= F && tab.back().first <= F))
		{
			auto work_tab = reverseTable(getNearest(reverseTable(shiftTable(tab, F)), 0, N));
			func.initNewtonCoef(work_tab, N);
		}
		vector<Functions::ld> ranges = Functions::separateRange(tab.front().second, tab.back().second, 1000);
		vector<Functions::ld> roots = func.roots(eps, ranges);
		for(auto elem: roots)
		{
			cout << "x = " << elem << endl << "|P(x) - F| = " << abs(func.evaluateNewtonInter(elem)) << endl;
		}
	}
	cout << endl;
	cout << endl;
}

void doRevInterpolation()
{
	cout << "Enter F" << endl << ">>";
	Functions::ld F = 0;
	cin >> F;
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
	Functions::ld eps = 0;
	cout << "Enter eps " <<  endl << ">>";
	cin >> eps;
	doFirstMethod(F, N);
	doSecondMethod(table[0].second, table[table.size() - 1].second, eps, N, F);
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
void printTable(vector<pair<Functions::ld, Functions::ld>> tab)
{
	for (int i = 0; i < tab.size(); ++i)
	{
		cout << "f(x_" << i << ") = " << tab[i].first << " | x_" << i << " = " << tab[i].second << endl;
	}
}

void printTable()
{
	for (int i = 0; i < table.size(); ++i)
	{
		cout << "f(x_" << i << ") = " << table[i].first << " | x_" << i << " = " << table[i].second << endl;
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
	cout << "The problem of reverse interpolation. Var-5. f(x) = 1 - exp(-2x)" << endl;
	cout << "1. Reverse Interpolation"<<endl;
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
			doRevInterpolation();
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
			break;
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
