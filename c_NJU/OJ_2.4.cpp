#include <iostream>
using namespace std;

int f(int n,int k);
int main()
{
	int n,k;
	while(cin >> n >> k)
	{
		cout <<f(n,k)<<endl;
	}
	return 0;
}

int f(int n,int k)
{
	int i;
	int a=0;
	if (n==0)
		return 1;
	else if (n==1)
		return 1;
	for (i=1;i<=k;i++)
	{
		if (n-i>=0)
			a += f(n-i,k);
	}	
	return a;
}
