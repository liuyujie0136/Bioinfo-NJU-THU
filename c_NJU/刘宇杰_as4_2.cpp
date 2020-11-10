#include <iostream>
#include <cmath> 
using namespace std;
int CountDigit(int n); //声明函数 
int KthDigit(int n, int k);//声明函数 

int main()
{
	int n,m;
	cout << "请输入正整数n:";
	cin >> n;
	m = CountDigit(n);//调用函数 
	cout <<"这是一个"<<m<<"位数。其右起第" << m/2 << "位上数字是:" << KthDigit(n,m/2) <<endl;//调用函数 
	return 0; 
}

int CountDigit(int n)//定义函数 
{
	int i;
	for (i=1;i<=n;i++)
	{
		n = n / 10;
		if (n < 10) break;
	}
	return i+1;
}

int KthDigit(int n, int k) //定义函数 
{
	n = n % (int)pow((double)10, k);
	n = n / pow((double)10, k-1);
	return n;
}


