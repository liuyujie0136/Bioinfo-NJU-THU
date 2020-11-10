#include <iostream>
using namespace std;

int compose(int n,int r)
{
	int i,a=1,b=1;
	for (i=r+1;i<=n;i++) a *= i;
	for (i=2;i<=n-r;i++) b *= i;
	return a/b;
}
int main()
{
	int n,r;
	cout << "本程序用于计算从n个数取r个数的所有选择个数。" << endl;
	cout <<"请输入n:" ;
	cin >> n;
	cout << "请输入r:";
	cin >> r;
	cout << "计算结果为:" << compose(n,r) <<endl;
	return 0;
}
