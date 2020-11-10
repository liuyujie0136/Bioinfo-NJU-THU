#include <iostream>
#include <cstring>
using namespace std;

int strlen(const char *s[]);
int strncmp(const char *s1[],const char *s2[],int n);
int myfind(char s1[],const char s2[]);
int main()
{
	char s1[100],s2[50];
	cin >> s1 >> s2;
	cout << myfind(s1,s2);
	return 0;
} 

int myfind(char s1[],const char s2[])
{
	int i,j;
	bool flag=true;
	int l1=strlen(s1),l2=strlen(s2);
	for (i=0;i<=(l1-l2);i++)
		if (strncmp(s1+i,s2,l2)==0)
		{
			flag=false;
			return i;
		}
	if (flag) return -1;
}
