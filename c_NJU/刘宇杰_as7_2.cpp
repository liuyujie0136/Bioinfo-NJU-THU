#include <iostream>
using namespace std;

int main()
{
	int a,b,c,m;
	cout << "请输入三个不同的正整数:" <<endl;
	cin >> a >> b >> c;
	if ((b<a&&a<c)||(c<a&&a<b))
		m = a;
	if ((a<b&&b<c)||(c<b&&b<a))
		m = b;
	if ((a<c&&c<b)||(b<c&&c<a))
		m = c;
	cout << "三个数中第二大的数是:" << m;
	return 0; 
}
