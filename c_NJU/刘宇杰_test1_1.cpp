#include <iostream>
using namespace std;

int main()
{
	int n,i,j,k,l;
	cout << "请输入一个正整数(0表示结束程序):";
	cin >> n;	
	while (n != 0)
	{
		for (i=0;i<=n;i++)
			for (j=0;j<=n;j++)
				for (k=0;k<=n;k++)
					for (l=0;l<=n;l++)
						if (n == i*i + j*j + k*k + l*l)
						{
							cout<<n<<"="<<i<<"*"<<i<<"+"<<j<<"*"<<j<<"+"<<k<<"*"<<k<<"+"<<l<<"*"<<l<<endl;
							goto End;
						}
		End:;			
		cout << "请输入一个正整数(0表示结束程序):";
		cin >> n;					
	}
	return 0;	
}
