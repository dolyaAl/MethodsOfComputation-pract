#include <iostream>
#include <iomanip>
#include "FunctionReader.h"
#include "Functions.h"

using namespace std;

void doCalc(long double A, long double B,long double eps, int N)
{
	Functions::Function a;
	cout << "roots count = " << a.countRoots(A, B, N) << endl;
	vector<Functions::ld> ranges = Functions::separateRange(A, B, N);
	vector<Functions::ld> roots = a.bisectionMethod(eps, ranges);
	cout << "bisection roots" << endl;
	cout << setprecision(15);
	for(auto elem: roots)
	{
		cout << elem << endl << "|f(x) - 0| = " << abs(a.evaluate(elem)) << endl;
	}
	cout << endl;
	roots = a.newtonMethod(eps, ranges);
	cout << "newton roots" << endl;
	for(auto elem: roots)
	{
		cout << elem << endl << "|f(x) - 0| = " << abs(a.evaluate(elem)) << endl;
	}
	cout << endl;
	roots = a.newtonModMethod(eps, ranges);
	cout << "newtonMod roots" << endl;
	for(auto elem: roots)
	{
		cout << elem << endl << "|f(x) - 0| = " << abs(a.evaluate(elem)) << endl;
	}
	cout << endl;
	roots = a.secantMethod(eps, ranges);
	cout << "secant roots" << endl;
	for(auto elem: roots)
	{
		cout << elem << endl << "|f(x) - 0| = " << abs(a.evaluate(elem)) << endl;
	}
	cout << endl;
}

int main()
{
	int N = 0;
	long double A = 0, B = 0, eps = 0;
	cout << "Enter N:" << endl << ">> ";
	cin >> N;
	cout << "Enter A:" << endl << ">> ";
	cin >> A;
	cout << "Enter B:" << endl << ">> ";
	cin >> B;
	cout << "Enter Eps:" << endl << ">> ";
	cin >> eps;
	doCalc(A, B, eps, N);
	return 0;
}
