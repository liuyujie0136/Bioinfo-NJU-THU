#include <iostream>
using namespace std;

int main()
{
	int a,b,c,m;
	cout << "������������ͬ��������:" <<endl;
	cin >> a >> b >> c;
	if ((b<a&&a<c)||(c<a&&a<b))
		m = a;
	if ((a<b&&b<c)||(c<b&&b<a))
		m = b;
	if ((a<c&&c<b)||(b<c&&c<a))
		m = c;
	cout << "�������еڶ��������:" << m;
	return 0; 
}
