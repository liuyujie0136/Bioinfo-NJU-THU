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
		/* ǰn-1��ȫװ�� ����ô����n����֮ǰ�κ�һ��������n����ȫװ��
		ǰn-1����һ��װ�ԣ���ô��n��ֻ�ܺ���һ�⽻�� ��
		ǰn-1����һ��װ�ԣ��൱��n-1(����λ����)���Ե�n-1��װ�ԣ��༴ǰn-2��ȫ��װ��*/ 
}
