#include <iostream>
using namespace std;

int main()
{
	int iYear,iMonth,iWeek,iDay,n,j;
	do
	{
		cout << "�������������:"; 
		cin >> iYear;
	}while (iYear%4 != 0 || iYear%100 == 0);
	do
	{
		cout << "�����������·�:"; 
		cin >> iMonth;
	}while (iMonth<1 || iMonth>12);
	
	cout << "\t\t" << iYear << " �� " << iMonth << " �� " << endl;
	cout << "��\tһ\t��\t��\t��\t��\t��" << endl;
	
	switch (iMonth)
	{
		case 1:n = 31;iYear--;iMonth = 13;break;
		case 3:case 5:case 7:case 8:case 10:case 12:n = 31;break;
		case 4:case 6:case 9:case 11:n = 30;break;
		case 2:n = 29;iYear--;iMonth = 14;break;
	}
	
	int c=iYear/100,y=iYear%100;

	iWeek = ((c/4)-2*c+y+(y/4)+(26*(iMonth+1)/10)) % 7; 
	for (j=1;j<iWeek+1;j++)
		cout << "\t";
	
	for (iDay=1;iDay<=n;iDay++)
	{
		iWeek = ((c/4)-2*c+y+(y/4)+(26*(iMonth+1)/10)+iDay-1) % 7; 
		cout << iDay <<"\t";
		if (iWeek==6)
			cout<<endl;
	}
	return 0;
}
