#include <iostream>
#include <cmath> 
using namespace std;
int CountDigit(int n); //�������� 
int KthDigit(int n, int k);//�������� 

int main()
{
	int n,m;
	cout << "������������n:";
	cin >> n;
	m = CountDigit(n);//���ú��� 
	cout <<"����һ��"<<m<<"λ�����������" << m/2 << "λ��������:" << KthDigit(n,m/2) <<endl;//���ú��� 
	return 0; 
}

int CountDigit(int n)//���庯�� 
{
	int i;
	for (i=1;i<=n;i++)
	{
		n = n / 10;
		if (n < 10) break;
	}
	return i+1;
}

int KthDigit(int n, int k) //���庯�� 
{
	n = n % (int)pow((double)10, k);
	n = n / pow((double)10, k-1);
	return n;
}


