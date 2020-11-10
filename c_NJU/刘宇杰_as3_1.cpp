#include <iostream>
using namespace std;

int main()
{
	int a,d,n,i,s = 0;
	cout << "请输入等差数列的项数n:";
	cin >> n;
	cout << "请输入等差数列的首项a:";
	cin >> a;
	cout << "请输入等差数列的公差d:";
	cin >> d;
	for (i = 1;i <= n;i++) 
	{
		s = s + a;
		a = a + d;
	}
	cout << "等差数列的和为:" << s << endl;
	return 0;
}
