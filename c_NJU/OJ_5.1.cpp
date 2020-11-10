#include <iostream>
#include <cstring> 
using namespace std;
#define N 257 

int strcmp(const char myin[],const char myout[]);

int main()
{
	int num=0,j;
	char myin[N]={""},myout[N]={""};
	cin >> myin;
	while (myin[num])
		num++;
	for (j=0;j<num;j++)
		myout[j]=myin[num-1-j];
	if (strcmp(myin,myout) == 0)
		cout <<"Y";
	else cout <<"N";
	return 0;
}
