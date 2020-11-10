#include <iostream>
using namespace std;

int IntRevs(int n);
int myPow(int n,int k);

int main()
{
	int n;
	cout << "请输入一个正整数:";
	cin >> n;
	cout << "该正整数的逆序数为:" << IntRevs(n) << endl;
	return 0;
}

int myPow(int n,int k)
{
	int i,p=1;
	for (i=1;i<=k;i++)
		p = p * n;
	return p;
}

int IntRevs(int n)
{
	int i,j,t,k=0,m=0;
	t = n;
	for (i=1;i<=t;i++)
	{
		t = t / 10;
		if (t < 10)  break;
	}
	k = i + 1;
	for (j=1;j<=i+1;j++)
	{
		m = m + ( n % 10) * myPow(10,k-1);
		n = n / 10;
		k = k - 1;  
	}
	return m;
}
