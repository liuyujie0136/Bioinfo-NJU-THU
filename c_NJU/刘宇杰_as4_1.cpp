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
	cout << "���������ڼ����n����ȡr����������ѡ�������" << endl;
	cout <<"������n:" ;
	cin >> n;
	cout << "������r:";
	cin >> r;
	cout << "������Ϊ:" << compose(n,r) <<endl;
	return 0;
}
