#include <iostream>
using namespace std;

int main()
{
	int n,i,j,k;
	while(cin >> n)
	{
		k=0;
		for (i=1;i<=n;i++)
		{
			j=i;
			while(j%5==0 && j!=0)
			{
				k++;
				j= j / 5;
			}
		}
		cout << k <<endl;
	}
	return 0;
}
