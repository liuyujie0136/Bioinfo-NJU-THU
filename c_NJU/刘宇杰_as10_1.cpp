#include <iostream>
using namespace std;
#define N 100

int main()
{
	int num=0,j;
	char myin[N],myout[N];
	cin >> myin;
	while (myin[num])
		num++;
	for (j=0;j<num;j++)
		myout[j]=myin[num-1-j];
	cout << myout;
	return 0;
}
