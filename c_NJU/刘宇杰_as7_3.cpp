#include <iostream>
#include <cmath>
using namespace std;

int main()
{
	int n,a,b,c,d,i;
	for (n=1100;n<=9988;n++)
	{
		a = n / 1000;
		b = n / 100 % 10;
		c = n % 100 / 10;
		d = n % 10;
		i = (int)sqrt(n);
		if (a==b && c==d && n==i*i)
		{
			cout << n;
			break;
		}
	}
	return 0;
}

