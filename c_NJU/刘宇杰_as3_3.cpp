#include <iostream>
#include <cmath> 
using namespace std;

int main()
{
	int n,i,j;
	cout << "请输入菱形行数:";
	cin >> n;
	for (i=1;i<=n;i++)
	{ 
	  for (j=1;j<=n;j++)
	  {
	  	if (j>=(n+1)/2-((n-1)/2-abs(i-(n+1)/2)) && j<=(n+1)/2+((n-1)/2-abs(i-(n+1)/2)))
	  	  cout << "*";
		else cout << " ";
	  }
	  cout << endl;
	}	
	return 0;
}
