#include <iostream>
#include <iomanip>
using namespace std;
int main()
{
	double a,b,c,d;
	cout<<"������һԪ���η���ax^2+b^x+c=0ϵ��a,b,c"<<endl;
	cin>>a>>b>>c;
	d=b*b-4*a*c;
	cout<<setiosflags(ios::fixed)<<setprecision(2)<<d<<endl;
	return 0;
}
