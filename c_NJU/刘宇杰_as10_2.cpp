#include <iostream>
#include <cstring>
using namespace std;
#define N 16

int strcmp(const char mykey[],const char key[]);

int main()
{
	int flag=0;
	char mykey[N],key[N]={"lyj987654"};
	cout << "Please input your key:";
	while (flag<3)
	{
		cin >> mykey;
		if (strcmp(mykey,key) == 0)
		{
			cout << "Pass!";
			break;
		}
		else if (flag <2)
			cout << "Wrong! Please try again:";
		flag++;
	}
	if (flag == 3)
		cout << "Your turns are used up. Good Luck!";
	return 0;
}
