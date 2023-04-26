#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include "Functions.h"

using namespace std;

using ld = long double;
Functions::Function func;

ld sincos(ld x)
{
	return sin(x);
}
ld Pervsincos(ld x)
{
	return -cos(x);
}

ld consta(ld x)
{
	return 5;
}
ld Pervconsta(ld x)
{
	return 5*x;
}

ld lin(ld x)
{
	return x*3/4;
}
ld Pervlin(ld x)
{
	return x*x*3/8;
}

ld quadr(ld x)
{
	return x*x+x;
}
ld Pervquadr(ld x)
{
	return x*x*x/3 + x*x/2;
}

ld threes(ld x)
{
	return x*x*x + x*x+x;
}
ld Pervthrees(ld x)
{
	return x*x*x*x/4 + x*x*x/3 + x*x/2;
}
ld e(ld x)
{
	return exp(2*x);
}
ld perve(ld x)
{
	return (1.0/2)*exp(2*x);
}
ld diffe(ld x, int n)
{
	return pow(2, n)*exp(2*x);
}


map<string, pair<ld(*)(ld), ld(*)(ld)>> evalf =
{
	{"sincos", {sincos, Pervsincos}},
	{"consta", {consta, Pervconsta}},
	{"lin", {lin, Pervlin}},
	{"quadr", {quadr, Pervquadr}},
	{"threes", {threes, Pervthrees}},
	{"e", {e, perve}}
};
map<string, string> nick_to_name =
{
	{"sincos", "sin(x) + cos(2x)"},
	{"consta", "5"},
	{"lin", "(3/4)x"},
	{"quadr", "x^2 + x"},
	{"threes", "x^3 + x^2 + x"},
	{"e", "e^x"}
};

bool exit(bool exit = false)
{
	static bool ex;
	if(exit)
	{
		ex = true;
	}
	return ex;
}

void do41(ld A, ld B, string func_s)
{
	func.setFunc(evalf[func_s].first);
	cout << "f(x) = " << nick_to_name[func_s] << "  A = " << A << "  B = " << B << endl;
	cout << "Left Rect" << endl;
	cout << "Numeric: " << func.quadrfleftRect(A, B) << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(func.quadrfleftRect(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	cout << "Right Rect" << endl;
	cout << "Numeric: " << func.quadrfrightRect(A, B) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.quadrfrightRect(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	cout << "Mid Rect" << endl;
	cout << "Numeric: " << func.quadrfmidRect(A, B) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.quadrfmidRect(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	cout << "Trapeze" << endl;
	cout << "Numeric: " << func.quadrfTrapeze(A, B) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.quadrfTrapeze(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	cout << "Simpson" << endl;
	cout << "Numeric: " << func.quadrfSimpson(A, B) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.quadrfSimpson(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	cout << "3/8" << endl;
	cout << "Numeric: " << func.quadrfThreeEighths(A, B) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.quadrfThreeEighths(A, B) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;

}
void do42(ld A, ld B, int M, string func_s)
{
	func.setFunc(evalf[func_s].first);
	ld theor = 0;
	int r = 0;
	cout << "f(x) = " << nick_to_name[func_s] << "  A = " << A << "  B = " << B << "  M = " << M << "  h = " << (B-A)/M << endl;
	
	cout << "Left Rect" <<  endl;
	cout << "Numeric: " << func.NIntleftRect(A, B, M) << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(func.NIntleftRect(A, B, M) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	if (func_s == "e")
	{
		r = 1;
		theor = (1.0 / 2) * diffe(B, r) * (B - A) * pow(((B - A) / M), r);
		cout << "Theor delta  =  " << theor << endl;
	}

	cout << "Right Rect" << endl;
	cout << "Numeric: " << func.NIntrightRect(A, B, M) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntrightRect(A, B, M) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	if (func_s == "e")
	{
		r = 1;
		theor = (1.0 / 2) * diffe(B, r) * (B - A) * pow(((B - A) / M), r);
		cout << "Theor delta  =  " << theor << endl;
	}

	cout << "Mid Rect" << endl;
	cout << "Numeric: " << func.NIntmidRect(A, B, M) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntmidRect(A, B, M) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	if (func_s == "e")
	{
		r = 2;
		theor = (1.0 / 24) * diffe(B, r) * (B - A) * pow(((B - A) / M), r);
		cout << "Theor delta  =  " << theor << endl;
	}

	cout << "Trapeze" << endl;
	cout << "Numeric: " << func.NIntTrapeze(A, B, M) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntTrapeze(A, B, M) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	if (func_s == "e")
	{
		r = 2;
		theor = (1.0 / 12) * diffe(B, r) * (B - A) * pow(((B - A) / M), r);
		cout << "Theor delta  =  " << theor << endl;
	}

	cout << "Simpson" << endl;
	cout << "Numeric: " << func.NIntSimpson(A, B, M) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntSimpson(A, B, M) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	if (func_s == "e")
	{
		r = 4;
		theor = (1.0 / 2880) * diffe(B, r) * (B - A) * pow(((B - A) / M), r);
		cout << "Theor delta  =  " << theor << endl;
	}
}
void do43(ld A, ld B, int M, ld L, string func_s)
{
	func.setFunc(evalf[func_s].first);
	do42(A, B, M, func_s);
	ld rung = 0;
	int r = 0;
	cout << endl << endl;
	cout << "f(x) = " << nick_to_name[func_s] << "  A = " << A << "  B = " << B << "  M*l = " << M*L << "  h = " << (B-A)/(M*L) << endl;

	cout << "Left Rect" << endl;
	cout << "Numeric: " << func.NIntleftRect(A, B, M*L) << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(func.NIntleftRect(A, B, M*L) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	r = 1;
	rung = (pow(L, r) * func.NIntleftRect(A, B, M * L) - func.NIntleftRect(A, B, M)) / (pow(L, r) - 1);
	cout << "Runge: " << rung << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(rung - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	
	cout << "Right Rect" << endl;
	cout << "Numeric: " << func.NIntrightRect(A, B, M*L) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntrightRect(A, B, M*L) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	r = 1;
	rung = (pow(L, r) * func.NIntrightRect(A, B, M * L) - func.NIntrightRect(A, B, M)) / (pow(L, r) - 1);
	cout << "Runge: " << rung << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(rung - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	

	cout << "Mid Rect" << endl;
	cout << "Numeric: " << func.NIntmidRect(A, B, M*L) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntmidRect(A, B, M*L) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	r = 2;
	rung = (pow(L, r) * func.NIntmidRect(A, B, M * L) - func.NIntmidRect(A, B, M)) / (pow(L, r) - 1);
	cout << "Runge: " << rung << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(rung - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	

	cout << "Trapeze" << endl;
	cout << "Numeric: " << func.NIntTrapeze(A, B, M*L) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntTrapeze(A, B, M*L) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	r = 2;
	rung = (pow(L, r) * func.NIntTrapeze(A, B, M * L) - func.NIntTrapeze(A, B, M)) / (pow(L, r) - 1);
	cout << "Runge: " << rung << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(rung - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	

	cout << "Simpson" << endl;
	cout << "Numeric: " << func.NIntSimpson(A, B, M*L) << "  Exact: " << (evalf[func_s].second(B) - evalf[func_s].second(A)) << "  delta = "  << abs(func.NIntSimpson(A, B, M*L) - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	r = 4;
	rung = (pow(L, r) * func.NIntSimpson(A, B, M * L) - func.NIntSimpson(A, B, M)) / (pow(L, r) - 1);
	cout << "Runge: " << rung << "  Exact: " << evalf[func_s].second(B) - evalf[func_s].second(A) << "  delta = "  << abs(rung - (evalf[func_s].second(B) - evalf[func_s].second(A))) << endl;
	
}

void printMainMenu()
{
	cout << "The problem of numeric integration." << endl;
	cout << "1. task 4.1"<<endl;
	cout << "2. task 4.2"<<endl;
	cout << "3. task 4.3"<<endl;
	cout << "4. Exit"<<endl;
}

void print41Menu() 
{
	cout << "Task 4.1" << endl;
	cout << "1. set A, B"<<endl;
	cout << "2. Eval integrals"<<endl;
	cout << "3. Change function"<<endl;
	cout << "4. Back"<<endl;
}
void print42Menu() 
{
	cout << "Task 4.2" << endl;
	cout << "1. set A, B"<<endl;
	cout << "2. set M"<<endl;
	cout << "3. Eval integrals"<<endl;
	cout << "4. Change function"<<endl;
	cout << "5. Back"<<endl;
}
void print43Menu() 
{
	cout << "Task 4.2" << endl;
	cout << "1. set A, B"<<endl;
	cout << "2. set M"<<endl;
	cout << "3. set l"<<endl;
	cout << "4. Eval integrals"<<endl;
	cout << "5. Change function"<<endl;
	cout << "6. Back"<<endl;
}
void printChangeFuncMenu() 
{
	cout << "Choose function" << endl;
	cout << "1. f(x) = sin(x)+cos(2x)"<<endl;
	cout << "2. f(x) = 5"<<endl;
	cout << "3. f(x) = (3/4)x"<<endl;
	cout << "4. f(x) = x^2 + x"<<endl;
	cout << "5. f(x) = x^3 + x^2 + x"<<endl;
	cout << "6. f(x) = e^x"<<endl;
	cout << "7. Back"<<endl;
}
string readChangecommand() 
{
	int com = 0;
	bool back = false;
	while (!back)
	{
		cin >> com;
		switch (com)
		{
		case 1:
		{
			return "sincos";
			break;
		}
		case 2:
		{
			return "consta";
			break;
		}
		case 3:
		{
			return "lin";
			break;
		}
		case 4:
		{
			return "quadr";
			break;
		}
		case 5:
		{
			return "threes";
			break;
		}
		case 6:
		{
			return "e";
			break;
		}
		case 7:
		{
			back = true;
			break;
		}
		default:
		{
			cout << "Try again" << endl;
		}
		}
	}
	return "";
}

void read41command() 
{
	int com = 0;
	bool back = false;
	ld A = 0;
	ld B = 0;
	string func_s = "sincos";
	while (!back)
	{
		system("cls");
		print41Menu();

		cin >> com;
		switch (com)
		{
		case 1:
		{
			cout << "A >> ";
			cin >> A;

			cout << "B >> ";
			cin >> B;
			break;
		}
		case 2:
		{
			do41(A, B, func_s);
			system("pause");
			break;
		}
		case 3:
		{
			system("cls");
			printChangeFuncMenu();
			func_s = readChangecommand();
			break;
		}
		case 4:
		{
			back = true;
			break;
		}
		default:
		{
			cout << "Try again" << endl;
		}
		}
	}

}
void read42command() 
{
	int com = 0;
	ld A = 0;
	ld B = 0;
	int M = 1;
	string func_s = "sincos";
	bool back = false;
	while (!back)
	{
		system("cls");
		print42Menu();
		cin >> com;
		switch (com)
		{
		case 1:
		{
			cout << "A >> ";
			cin >> A;

			cout << "B >> ";
			cin >> B;
			break;
		}
		case 2:
		{
			cout << "M >> ";
			cin >> M;
			break;
		}
		case 3:
		{
			do42(A, B, M, func_s);
			system("pause");
			break;
		}
		case 4:
		{
			system("cls");
			printChangeFuncMenu();
			func_s = readChangecommand();
			break;
		}
		case 5:
		{
			back = true;
			break;
		}
		default:
		{
			cout << "Try again" << endl;
		}
		}
	}
}
void read43command() 
{
	int com = 0;
	ld A = 0;
	ld B = 0;
	int M = 1;
	ld L = 1;
	string func_s = "sincos";
	bool back = false;
	while (!back)
	{
		system("cls");
		print43Menu();
		cin >> com;
		switch (com)
		{
		case 1:
		{
			cout << "A >> ";
			cin >> A;

			cout << "B >> ";
			cin >> B;
			break;
		}
		case 2:
		{
			cout << "M >> ";
			cin >> M;
			break;
		}
		case 3:
		{
			cout << "L >> ";
			cin >> L;
			break;
		}
		case 4:
		{
			do43(A, B, M, L, func_s);
			system("pause");
			break;
		}
		case 5:
		{
			system("cls");
			printChangeFuncMenu();
			func_s = readChangecommand();
			break;
		}
		case 6:
		{
			back = true;
			break;
		}
		default:
		{
			cout << "Try again" << endl;
		}
		}
	}
}

void readMainCommand()
{
	int N = 0;
	cin >> N;
	switch(N)
	{
		case 1:
		{
			read41command();
			system("pause");
			break;
		}
		case 2:
		{
			read42command();
			break;
		}
		case 3:
		{
			read43command();
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
		printMainMenu();
		readMainCommand();
		system("cls");
	}
	return 0;
}
