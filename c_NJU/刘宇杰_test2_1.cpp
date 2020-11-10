#include <iostream>
using namespace std;

double e(int x);
int main()
{
	int x;
	cin >> x;
	cout << e(x);
	return 0;
}

double e(int x)
{
	if (x == 0) 
		return 1.0;
	else
	{
		int i=1,j;
		double sum=1.0,temp;
		do
		{
			temp = 1.0;
			for (j=1;j<=i;j++)
				temp = temp*double(x)/double(j);
			sum += temp;
			i += 1;
		}while (temp>=1e-7);
		return sum;
	}
}
