#include <iostream>
using namespace std;

int main()
{
	int n,i,a1,a2,t;
	cout << "������һ��̨����(��������0��ʾ��������):";
	cin >> n;	
	while (n !=0)
	{
		a1=1;a2=1;
		for (i=2;i<=n;i++)
		{
			t=a1+a2;
			a1=a2;
			a2=t;
		}
		cout << "��������Ϊ:" << a2 << endl << "������һ��̨����(��������0��ʾ��������):";	
		cin >> n;					
	}
	return 0;		
}

