#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

int main()
{
	double x0,x1,x2=0;
	cout << "������һ����:";
	cin >> x0;
	x1 = x0;
	x2 = (2*x1+x0/pow(x1,2))/3;
	while (abs(x2-x1) > 0.000001)
	{
		x1 = x2;
		x2 = (2*x1+x0/pow(x1,2))/3;
	}
	cout << x0 << "��������Ϊ:" << fixed << setprecision(2) << x2;
	return 0; 
}
