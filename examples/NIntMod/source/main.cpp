#include <iostream>
#include <iomanip>
#include "Functions.h"

using namespace std;
using ld = long double;
void do51()
{
	string fu_str = "";
	string we_str = "";
	ld A = 0;
	ld B = 0;
	ld tochn = 0;
	int n = 0;
	cout << "f(x) = ";
	cin >> fu_str;
	cout << "ro(x) = ";
	cin >> we_str;
	cout << "N = ";
	cin >> n;
	cout << "A = ";
	cin >> A;
	cout << "B = ";
	cin >> B;
	Functions::Function fu2(we_str + "*" + fu_str);
	for (int i = 0; i < 10000; ++i)
	{
		tochn+= fu2.quadrfSimpson(i * ((B-A) / 10000), (i + 1) * ((B-A) / 10000));
	}
	Functions::Function fu(fu_str);
	Functions::Weight we(we_str);
	auto moms = we.get_weight_moms(A, B, 2 * n);
	for (int i = 0; i < 2*n; ++i)
	{
		cout << "mu_" << to_string(i) << " = " << moms[i] << endl;
	}
	cout << endl;
	ld integr = fu.NIntInterQF(A, B, we, n);
	cout << "NIntINT[" << fu_str << ", " << A << ", " << B << "] = " << integr << endl << "delta = " << abs(tochn - integr);
	cout << endl;
	integr = fu.NIntMAGQF(A, B, we, n);
	cout << "NIntMAG[" << fu_str << ", " << A << ", " << B << "] = " << integr << endl << "delta = " << abs(tochn - integr);
	cout << endl;
}
void do52gauss()
{
	string fu_str = "";
	string we_str = "";
	ld A = 0;
	ld B = 0;
	ld tochn = 0;
	int n1g	= 0;
	int n2g	= 0;
	int n3g	= 0;
	int n4g	= 0;
	cout << "f(x) = ";
	cin >> fu_str;
	cout << " A = ";
	cin >> A;
	cout << " B = ";
	cin >> B;
	cout << " N1 = ";
	cin >> n1g;
	cout << " N2 = ";
	cin >> n2g;
	cout << " N3 = ";
	cin >> n3g;
	cout << " N4 = ";
	cin >> n4g;
	cout << " tochn = ";
	cin >> tochn;
	Functions::Function fu(fu_str);
	auto gau = fu.NIntGauss(A, B, n1g, true);
	cout << endl << "NInt[" << fu_str << ", " << A << ", " << B << ", " << n1g << "] = " << gau << endl << "delta = " << abs(tochn - gau) << endl << endl;
	gau = fu.NIntGauss(A, B, n2g, true);
	cout << endl << "NInt[" << fu_str << ", " << A << ", " << B << ", " << n2g << "] = " << gau << endl << "delta = " << abs(tochn - gau) << endl << endl;
	gau = fu.NIntGauss(A, B, n3g, true);
	cout << endl << "NInt[" << fu_str << ", " << A << ", " << B << ", " << n3g << "] = " << gau << endl << "delta = " << abs(tochn - gau) << endl << endl;
	gau = fu.NIntGauss(A, B, n4g, true);
	cout << endl << "NInt[" << fu_str << ", " << A << ", " << B << ", " << n4g << "] = " << gau << endl << "delta = " << abs(tochn - gau) << endl << endl;
	cout <<  endl;
}
void do52meler()
{
	string fu_str = "";
	string we_str = "";
	int n1m = 0;
	int n2m = 0;
	int n3m = 0;
	int n4m = 0;
	ld tochn = 0;
	cout << "f(x) = ";
	cin >> fu_str;
	cout << "N1 = ";
	cin >> n1m;
	cout << "N2 = ";
	cin >> n2m;
	cout << "N3 = ";
	cin >> n3m;
	cout << "N4 = ";
	cin >> n4m;
	Functions::Function fu(fu_str);
	auto meler = fu.NIntMehler(n1m);
	cout << "NInt[" << fu_str << ", " << n1m << "] = " << meler << endl << "delta = " << abs(tochn - meler) << endl << endl;
	meler = fu.NIntMehler(n2m);
	cout << "NInt[" << fu_str << ", " << n2m << "] = " << meler << endl << "delta = " << abs(tochn - meler) << endl << endl;
	meler = fu.NIntMehler(n3m);
	cout << "NInt[" << fu_str << ", " << n3m << "] = " << meler << endl << "delta = " << abs(tochn - meler) << endl << endl;
	meler = fu.NIntMehler(n4m);
	cout << "NInt[" << fu_str << ", " << n4m << "] = " << meler << endl << "delta = " << abs(tochn - meler) << endl << endl;
	cout <<  endl;

}
void do53()
{
	string fu_str = "";
	string we_str = "";
	ld A = 0;
	ld B = 0;
	ld tochn = 0;
	int n = 0;
	int m = 0;
	cout << "f(x) = ";
	cin >> fu_str;
	cout << "A = ";
	cin >> A;
	cout << "B = ";
	cin >> B;
	cout << "N = ";
	cin >> n;
	cout << "M = ";
	cin >> m;
	cout << "tochn = ";
	cin >> tochn;
	Functions::Function fu(fu_str);
	ld ans = 0;
	for (int i = 0; i < m; ++i)
	{
		ans += fu.NIntGauss(A + i*(B-A)/m , A+(i+1)*(B-A)/m, n, i==0);
	}
	cout << endl << "NInt[f(x), A, B] = " << ans << endl << "delta = " << abs(tochn - ans) << endl << endl;
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
	cout << "The problem of numerical integration." << endl;
	cout << "1. 5.1"<<endl;
	cout << "2. 5.2 GAUSS"<<endl;
	cout << "3. 5.2 MELER"<<endl;
	cout << "4. 5.3"<<endl;
	cout << "5. Exit"<<endl;
}

void readCommand()
{
	int N = 0;
	cin >> N;
	switch(N)
	{
		case 1:
		{
			system("cls");
			do51();
			system("pause");
			break;
		}
		case 2:
		{
			system("cls");
			do52gauss();
			system("pause");
			break;
		}
		case 3:
		{
			system("cls");
			do52meler();
			system("pause");
			break;
		}
		case 4:
		{
			system("cls");
			do53();
			system("pause");
			break;
		}
		case 5:
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
	/*Functions::Function fu("x^1");
	Functions::Weight weight("1");
	int N = 0;
	while (true)
	{
		cin >> N;
		for (int i = 1; i < 2*N; ++i)
		{
			string a = "x^" + to_string(i);
			fu.updateFunc("x");
			cout <<  fu.NIntMehler(N) << endl;
		}
	}
	*/
	cout << setprecision(15);
	while(!exit())
	{
		printMenu();
		readCommand();
		system("cls");
	}
	return 0;
}
