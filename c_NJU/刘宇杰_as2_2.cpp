#include <iostream>
#include <iomanip>
using namespace std;
int main()
{
	double x,y;
	cout << "�����빫����x" << endl;
	cin >> x;
	if (x <= 3)
	  y = 9 + 1;
	else
	  y = 9 + 1 + 2.4 * (x-3);
	cout << "Ӧ����" << fixed << setprecision(1) << y << "Ԫ" << endl;
	return 0; 

}
