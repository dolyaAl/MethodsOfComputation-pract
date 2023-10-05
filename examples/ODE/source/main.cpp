#include <iostream>
#include <iomanip>
#include <cmath>
#include "Functions.h"

using namespace std;
using ld = long double;
vector<ld> tochn;
vector<ld> Teilor;
vector<ld> Adams;
vector<ld> Runge;
vector<ld> Euler;
vector<ld> Euler1;
vector<ld> Euler2;
vector<vector<ld>> qm;
void initToch(string ans, int n, ld h, ld x0)
{
	FunctionParser_ld pars;
	pars.AddConstant("PI", acos(-1));
	pars.Parse(ans, "x");
	long double x = 0;
	for (int i = -2; i < n+1; ++i)
	{
		x = x0 + h * i;
		tochn.push_back(pars.Eval(&x));
	}
}
void printToch()
{
	for (int i = 0; i < tochn.size(); ++i)
	{
		cout << "x0 + h*(" << to_string(i - 2) << ") = " << tochn[i] << endl;
	}
}
vector<string> getDfs(string func, int n)
{
	static int N = -1;
	static vector<string> dfs = {"y", func};
	if (N == n)
	{
		return dfs;
	}
	else
	{
		N = n;
	}
	Ev3::ExpressionParser pars;
	pars.SetVariableID("x", 0);
	pars.SetVariableID("y", 1);
	int nerr = 0;

	Ev3::Expression exp1;
	string gx = "";
	string gy = "";
	for (int i = dfs.size(); i <= n; ++i)
	{
		exp1 = pars.Parse(dfs.back().c_str(), nerr);
		gx = Ev3::Diff(exp1, 0)->ToString();
		gy = Ev3::Diff(exp1, 1)->ToString();
		dfs.push_back(gx + "+(" + gy + ")*(" + func + ")");
	}
	return dfs;
}
string doTeilor(string func, ld x0, ld y0)
{
	auto dfs = getDfs(func, 10);
	FunctionParser_ld parser;
	stringstream stream;
	stream << setprecision(15);
	string y = "";
	ld vars[2] = { x0, y0 };
	for (int i = 0; i < dfs.size(); ++i)
	{
		parser.Parse(dfs[i], "x,y");
		stream << "(" << parser.Eval(vars) << ")" << "*((x - " << x0 << ")^(" << i << "))/(" << Functions::factorial(i) << ")+";
	}
	y = stream.str();
	y.pop_back();
	return y;
}

void initQmByTeil(string func, ld h, ld x0)
{
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	qm.push_back(vector<ld>());
	for (int i = -2; i < 3; ++i)
	{
		ld var[2] = {x0 + h * i, Teilor[i + 2]};
		qm[0].push_back(pars.Eval(var));
	}
	for(int i = 1; i < 5; ++i)
	{
		qm.push_back(vector<ld>());
		for (int j = 0; j+i < 5; ++j)
		{
			qm[i].push_back(qm[i - 1][j + 1] - qm[i - 1][j]);
		}
	}
}
void doAdams(string func, ld h, ld x0,int n)
{
	initQmByTeil(func, h, x0);
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	for (int i = -2; i < 3; ++i)
	{
		Adams.push_back(Teilor[i + 2]);
	}
	ld tmp = 0;
	ld ym1 = 0;
	for (int i = 2; i < n; ++i)
	{
		tmp = qm[0][i + 2] + qm[1][i + 1] / 2 + 5 * qm[2][i] / 12 + 3 * qm[3][i - 1] / 8 + 251 * qm[4][i - 2] / 720;
		ym1 = Adams[i+2] + h*tmp;
		Adams.push_back(ym1);
		ld var[2] = {x0 + h * (i+1), ym1};
		qm[0].push_back(pars.Eval(var));
		for(int j = 1; j < 5; ++j)
		{
			ld qm1 = qm[j - 1][qm[j - 1].size() - 1];
			ld qm2 = qm[j - 1][qm[j - 1].size() - 2];
			qm[j].push_back(qm1 - qm2);
		}
	}
}
void doRunge(string func, ld h, ld x0, ld y0, int n)
{
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	Runge.push_back(y0);
	ld k1 = 0;
	ld k2 = 0;
	ld k3 = 0;
	ld k4 = 0;
	ld vars[2];
	for (int i = 0; i < n; ++i)
	{
		vars[0] = x0 + h*i; vars[1] = Runge.back();
		k1 = h * pars.Eval(vars);
		vars[0] = x0 + h*i + h/2; vars[1] = Runge.back() + k1/2;
		k2 = h * pars.Eval(vars);
		vars[0] = x0 + h*i + h/2; vars[1] = Runge.back() + k2/2;
		k3 = h * pars.Eval(vars);
		vars[0] = x0 + h*i + h; vars[1] = Runge.back() + k3;
		k4 = h * pars.Eval(vars);
		Runge.push_back(Runge.back() + (k1 + 2 * k2 + 2 * k3 + k4) / 6);
	}

}
void doEuler(string func, ld h, ld x0, ld y0, int n)
{
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	Euler.push_back(y0);
	ld vars[2];
	for (int i = 0; i < n; ++i)
	{
		vars[0] = x0 + h * i; vars[1] = Euler.back();
		Euler.push_back(Euler.back() + h*pars.Eval(vars));
	}
}
void doEuler1(string func, ld h, ld x0, ld y0, int n)
{
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	Euler1.push_back(y0);
	ld vars[2];
	for (int i = 0; i < n; ++i)
	{
		vars[0] = x0 + h * i; vars[1] = Euler1.back();
		vars[1] = Euler1.back() + h * pars.Eval(vars) / 2;
		vars[0] = x0 + h * i + h/2; 
		Euler1.push_back(Euler1.back() + h*pars.Eval(vars));
	}

}
void doEuler2(string func, ld h, ld x0, ld y0, int n)
{
	FunctionParser_ld pars;
	pars.Parse(func, "x,y");
	Euler2.push_back(y0);
	ld vars1[2];
	ld vars2[2];
	for (int i = 0; i < n; ++i)
	{
		vars1[0] = x0 + h * i; 
		vars1[1] = Euler2.back();
		vars2[0] = x0 + h * (i + 1); 
		vars2[1] = Euler2.back() + h * pars.Eval(vars1);
		Euler2.push_back(Euler2.back() + h*(pars.Eval(vars1)+pars.Eval(vars2))/2);
	}
}
void readData(string& fu_str, string& ans_str, int& n, ld& h, ld& y0, ld& x0)
{
	cout << "f(x, y) = ";
	cin >> fu_str;
	cout << "x0 = ";
	cin >> x0;
	cout << "y(x0) = ";
	cin >> y0;
	cout << "y(x) = ";
	cin >> ans_str;
	cout << "N = ";
	cin >> n;
	cout << "h = ";
	cin >> h;
}
void cl()
{
	tochn.clear();
	Teilor.clear();
	Adams.clear();
	qm.clear();
	Runge.clear();
	Euler.clear();
	Euler1.clear();
	Euler2.clear();
}
void doit()
{
	cl();
	string fu_str = "";
	string ans_str= "";
	int n = 0;
	ld h = 0;
	ld y0 = 0;
	ld x0 = 0;
	readData(fu_str, ans_str, n, h, y0, x0);
	initToch(ans_str, n, h, x0);
//	printToch();
	auto T = doTeilor(fu_str, x0, y0);
	FunctionParser_ld parser;
	parser.Parse(T, "x");
	ld x = 0;
	for (int i = -2; i < n+1; ++i)
	{
		x = x0 + h * i;
		Teilor.push_back(parser.Eval(&x));
	}
	/*cout << "######### TEILOR #########" << endl;
	for (int i = 0; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " << abs(tochn[i] - Teilor[i]) << endl;
	}
	cout << endl;*/
	doAdams(fu_str, h, x0, n);
	/*cout << "######### ADAMS #########" << endl;
	for (int i = 0; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " <<abs(tochn[i] - Adams[i]) << endl;
	}
	cout << endl;*/
	doRunge(fu_str, h, x0, y0, n);
	/*cout << "######### RUNGE #########" << endl;
	for (int i = 2; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " <<abs(tochn[i] - Runge[i-2]) << endl;
	}
	cout << endl;*/
	doEuler(fu_str, h, x0, y0, n);
	/*cout << "######### EULER #########" << endl;
	for (int i = 2; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " <<abs(tochn[i] - Euler[i-2]) << endl;
	}
	cout << endl;*/
	doEuler1(fu_str, h, x0, y0, n);
	/*cout << "######### EULER1 #########" << endl;
	for (int i = 2; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " <<abs(tochn[i] - Euler1[i-2]) << endl;
	}
	cout << endl;*/
	doEuler2(fu_str, h, x0, y0, n);
	/*cout << "######### EULER2 #########" << endl;
	for (int i = 2; i < n + 3; ++i)
	{
		cout << "dy(" << x0 + h * (i - 2) << ") = " <<abs(tochn[i] - Euler2[i-2]) << endl;
	}
	cout << endl;*/
	cout << "TE dyN = " <<abs(tochn[n+2] - Teilor[n+2]) << endl;
	cout << "AD dyN = " <<abs(tochn[n+2] - Adams[n+2]) << endl;
	cout << "RU dyN = " <<abs(tochn[n+2] - Runge[n]) << endl;
	cout << "E0 dyN = " <<abs(tochn[n+2] - Euler[n]) << endl;
	cout << "E1 dyN = " <<abs(tochn[n+2] - Euler1[n]) << endl;
	cout << "E2 dyN = " <<abs(tochn[n+2] - Euler2[n]) << endl;
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
	cout << "The problem of numerical solve DE." << endl;
	cout << "1. solve problem"<<endl;
	cout << "2. Exit"<<endl;
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
			doit();
			system("pause");
			break;
		}
		case 2:
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
