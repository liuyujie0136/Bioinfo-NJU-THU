#include <iostream>
using namespace std;

bool zhishu(int i);

int main()
{
	int a,i,j;
	cin >> a;
	for (i=2;i<=a;i++) 
	{
		j = 0;
		if (zhishu(i))
		{
			while (a % i == 0)
			{
				j = j + 1;
				a = a / i;
			}
			if (j != 0)
				cout << i << "(" << j << ")";
		}
	}
	return 0;
}

bool zhishu(int i)
{
	int k;
	bool flag=true;
	for (k=2;k<i;k++)
		if (i % k == 0)
		{
			flag = false;
			break;
		}	
	return flag;
}
