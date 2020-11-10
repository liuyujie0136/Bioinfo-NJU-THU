#include <iostream>
using namespace std;

int main()
{
	double i,s=0.0,t=-1.0;
	for (i=1.0;i<=100.0;i++)
	{
		t=-t;
		s+=t*1.0/i;
	}
	cout<<s;
	return 0;
}
