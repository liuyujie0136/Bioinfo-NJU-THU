#include <iostream>
using namespace std;
#define N 100000

int mysort(int *p,int *q);
int main()
{
	int n,i,j,*p,*q,a[N];
	cout << "Please imput the number of int. :";
	cin >> n;
	cout <<"Please imput the numbers:"<<endl;
	for (i=0;i<n;i++)
		cin >> a[i];
	for (i=0;i<n;i++)
	{
		p=&a[i];
		for (j=i+1;j<n;j++)
		{
			q=&a[j];
			mysort(p,q);
		}
	}
	for (i=0;i<n;i++)
		cout << a[i] <<" ";
	return 0;
}

int mysort(int *p,int *q)
{
	int t;
	if (*p > *q)
	{
		t=*p;
		*p=*q;
		*q=t;
	}
}
