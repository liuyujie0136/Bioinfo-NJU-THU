#include <iostream>
using namespace std;
#define RATIO 4
int main()
{
	char aaa[]="*",aab[]=" ";
	int n,i,j,k,a;
	cin >> n;
	for (i=0;i<=2*n;i++)
	{
		for (j=0-n;j<=n;j++)
			if ((RATIO*(i-n)*(i-n)+(j-n)*(j-n))<=RATIO*n*n)
			 {
			     cout << aaa;
			     a = j;
			     break;
			 } 
			else cout <<aab;
		for (k=j;k<=n;k++)
			cout <<aab;
		for (j=n;j<=4*n;j++)
			if (2*n-j==a)
			 {
			     cout <<aaa;
			 } 
			else cout <<aab;
		cout << endl;
	}
	return 0;
}
