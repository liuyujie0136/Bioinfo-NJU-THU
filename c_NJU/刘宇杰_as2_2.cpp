#include <iostream>
#include <iomanip>
using namespace std;
int main()
{
	double x,y;
	cout << "请输入公里数x" << endl;
	cin >> x;
	if (x <= 3)
	  y = 9 + 1;
	else
	  y = 9 + 1 + 2.4 * (x-3);
	cout << "应付款" << fixed << setprecision(1) << y << "元" << endl;
	return 0; 

}
