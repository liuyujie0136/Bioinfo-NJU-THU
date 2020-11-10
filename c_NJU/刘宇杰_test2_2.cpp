#include <iostream>
using namespace std;

int count(int n);
int main()
{
	int n;
	cin >> n;
	cout << count(n);
	return 0;
}

int count(int n)
{
	int i,j,k,s=1,t=1;
	if (n == 1)
		return 0;
	else if (n == 2)
		return 1;
	else
		return (n-1)*(count(n-1)+count(n-2));
		/* 前n-1封全装错 ，那么将第n封与之前任何一个交换，n个将全装错；
		前n-1个有一个装对，那么第n封只能和这一封交换 ；
		前n-1个有一个装对，相当于n-1(空闲位置数)乘以第n-1个装对，亦即前n-2个全部装错；*/ 
}
