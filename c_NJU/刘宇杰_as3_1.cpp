#include <iostream>
using namespace std;

int main()
{
	int a,d,n,i,s = 0;
	cout << "������Ȳ����е�����n:";
	cin >> n;
	cout << "������Ȳ����е�����a:";
	cin >> a;
	cout << "������Ȳ����еĹ���d:";
	cin >> d;
	for (i = 1;i <= n;i++) 
	{
		s = s + a;
		a = a + d;
	}
	cout << "�Ȳ����еĺ�Ϊ:" << s << endl;
	return 0;
}
