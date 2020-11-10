#include <iostream>
using namespace std;

int IntRevsR(int n);
int myPow(int n,int k);

int main()
{
	int n;
	cout << "请输入一个正整数:";
	cin >> n;
	cout << "该正整数的逆序数为:" << IntRevsR(n) << endl;
	return 0;
}

int myPow(int n,int k)
{
	int i,p=1;
	for (i=1;i<=k;i++)
		p = p * n;
	return p;
}

int IntRevsR(int n)
{
	int i,k,m,t=n;
	for (i=1;i<=t;i++)
	{
		t = t / 10;
		if (t < 10)  break;
	}
	k = i + 1;
	if ( n < 10 )
		return n;
	else
	{
		m = n % 10;
		n = n / 10;	
		return IntRevsR(n)+m*myPow(10,k-1);
	}
}
