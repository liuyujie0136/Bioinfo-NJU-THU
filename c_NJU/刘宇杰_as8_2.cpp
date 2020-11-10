#include <iostream>
using namespace std;
#define N 100

int main()
{
	int i=0,j,c=0,l=0,d=0;
	char a[N];
	int location[N][3];
	while(cin >> a[i])
	{
		if (a[i]>='A' && a[i]<='Z')
		{
			location[c][0] = i;
			c += 1;
		}
		
		if (a[i]>='a' && a[i]<='z')
		{
			location[l][1] = i;
			l += 1;
		}
		
		if (a[i]>='0' && a[i]<='9')
		{
			location[d][2] = i;
			d += 1;
		}
		i++;
	}
	
	cout << "Capital num: " << c <<": ";
	for (j=0;j<c;j++)
		cout << a[location[j][0]] <<" ";
	cout << endl;
	
	cout << "Lowercase num: " << l <<": ";
	for (j=0;j<l;j++)
		cout << a[location[j][1]] <<" ";
	cout << endl;
	
	cout << "Digits num: " << d <<": ";
	for (j=0;j<d;j++)
		cout << a[location[j][2]] <<" ";
	
	return 0;
}
