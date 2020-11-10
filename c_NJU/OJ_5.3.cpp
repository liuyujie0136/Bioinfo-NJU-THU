#include <iostream>
using namespace std;
#define N 10000

int main()
{
	char key[N][60]={""};
	int i,a,b,c,d,n=0;
	cin >> key[n];
	while (cin)
	{
		i=0;a=0;b=0;c=0;d=0;
		while (key[n][i])
		{
			if (key[n][i]>='A' && key[n][i]<='Z')
				a=1;
			else if (key[n][i]>='a' && key[n][i]<='z')
				b=1;
			else if (key[n][i]>='0' && key[n][i]<='9')
				c=1;
			else if (key[n][i]=='~'||key[n][i]=='!'||key[n][i]=='@'||key[n][i]=='#'||key[n][i]=='$'||key[n][i]=='%'||key[n][i]=='^')
				d=1;
			i++;
		}
		if ((a+b+c+d)>=3 && i>=9 && i<=17)
			cout << "YES"<<endl;
		else cout << "NO"<<endl;
		n++;
		cin >> key[n];
	}
	return 0;
}
