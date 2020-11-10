#include <iostream>
using namespace std;

int main()
{
	int n,i,a1,a2,t;
	cout << "请输入一个台阶数(正整数，0表示结束程序):";
	cin >> n;	
	while (n !=0)
	{
		a1=1;a2=1;
		for (i=2;i<=n;i++)
		{
			t=a1+a2;
			a1=a2;
			a2=t;
		}
		cout << "跳法总数为:" << a2 << endl << "请输入一个台阶数(正整数，0表示结束程序):";	
		cin >> n;					
	}
	return 0;		
}

