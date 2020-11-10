#include <iostream>
using namespace std;

int main()
{
	int *p=new int [10],n,i,count=0,max=10;
	cin >> n;
	while (n!=-1) //using -1 to stop the imput
	{
		p[count]=n;
		count++;
		if (count>=max)
		{
			max+=10;
			int *q=new int[max];
			for (i=0;i<count;i++)
				q[i]=p[i];
			delete []p;
			p = q;
		}
		cin >>n;
	}
	
	for (i=0;i<count;i++)
		cout << p[i] << " ";
	delete []p;
	return 0;
}
