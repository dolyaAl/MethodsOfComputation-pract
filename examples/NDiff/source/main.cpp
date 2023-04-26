#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "Functions.h"

using namespace std;

vector<pair<Functions::ld, Functions::ld>> table;
vector<vector<Functions::ld>> diff_tab;

Functions::ld f(Functions::ld x)
{
	return exp(1.5*x);
}

Functions::ld df(Functions::ld x)
{
	return 1.5 * f(x);
}

Functions::ld d2f(Functions::ld x)
{
	return 1.5 * 1.5 * f(x);
}

void doNumDiff()
{
	if(!diff_tab.empty())
	{
		diff_tab.clear();
	}
	using ld = Functions::ld;
	ld h = table[1].second - table[0].second;
	ld x = 0;
	ld fun = 0;
	ld dfun = 0;
	ld d2fun = 0;
	for (int i = 0; i < table.size(); ++i)
	{
		x = table[i].second;
		fun = table[i].first;
		if(i == 0)
		{
			dfun = (-3 * table[0].first + 4 * table[1].first - table[2].first) / (2 * h);
			d2fun = (2 * table[0].first - 5 * table[1].first + 4 * table[2].first - table[3].first) / (h * h);
		}
		else if(i == table.size() - 1)
		{
			dfun = (3 * table[i].first - 4 * table[i - 1].first + table[i-2].first) / (2 * h);
			d2fun = (2 * table[i].first - 5 * table[i - 1].first + 4 * table[i - 2].first - table[i - 3].first) / (h * h);
		}
		else
		{
			dfun = (table[i+1].first - table[i - 1].first) / (2 * h);
			d2fun = (table[i+1].first - 2*table[i].first + table[i-1].first) / (h * h);
		}
		diff_tab.push_back({ x, fun, dfun, abs(dfun - df(x)), d2fun, abs(d2fun - d2f(x)) });
	}
}

void initTable()
{
	if(!table.empty())
	{
		table.clear();
	}
	Functions::ld A = 0;
	Functions::ld h = 0;
	Functions::ld x = 0;
	int m = 0;
	cout << "Enter M (units count)" << endl << ">>";
	cin >> m;
	cout << "Enter left edge" << endl << ">>";
	cin >> A;
	cout << "Enter h edge" << endl << ">>";
	cin >> h;
	--m;
	for (int i = 0; i <= m; ++i)
	{
		x = A + h*i;
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

void printDiffTable()
{
	ofstream fout("difftab.txt");
	for (int i = 0; i < table.size(); ++i)
	{
		fout  <<  " x_" << i << " = " << diff_tab[i][0] << "\t***\tf(x_" << i << ") = "<< diff_tab[i][1]<<
			" \t***\tf'num(x_" << i << ") = "<< diff_tab[i][2]<<" \t***\t|f'num(x_" << i << ") - f'teor(x_" << i << ")| = "<< diff_tab[i][3]<<
			" \t***\tf''num(x_" << i << ") = "<< diff_tab[i][4]<<" \t***\t|f''num(x_" << i << ") - f''teor(x_" << i << ")| = "<< diff_tab[i][5]<<endl;
	}
	fout.close();
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
	cout << "The problem of Numeric differentiation. Var-5. f(x) = exp(3/2 * x)" << endl;
	cout << "1. Numeric differentiation"<<endl;
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
			doNumDiff();
			printDiffTable();
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
