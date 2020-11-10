#include <iostream>
#include <cmath>
using namespace std;
int main()
{
	int n,m;
	cout << "Please input n" << endl;
	cin >> n;
	m = abs(1-n);
    if (m>n)
	  cout << "|1-" << n << "| > " << n << "-1" << endl;
	else
	  cout << "|1-" << n << "| = " << n << "-1" << endl;
	return 0;
}
