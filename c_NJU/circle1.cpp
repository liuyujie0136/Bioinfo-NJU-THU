#include <iostream>
using namespace std;

int main()
{
	int n,i,j,k,a;
	cin >> n;
	for (i=0;i<=2*n;i++)
	{
		for (j=0-n;j<=n;j++)
			if ((4*(i-n)*(i-n)+(j-n)*(j-n))<=4*n*n)
			 {
			     cout << "*";
			     a = j;
			     break;
			 } 
			else cout <<" ";
		for (k=j;k<=n;k++)
			cout <<" ";
		for (j=n;j<=3*n;j++)
			if (2*n-j==a)
			 {
			     cout <<"*";
			 } 
			else cout <<" ";
		cout << endl;
	}
	return 0;
}
