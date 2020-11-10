#include <iostream>
using namespace std;
#define N 100

int main()
{
	int i=0,j=0,t,p,r,imax;
	double s=0.0,max;
	int flag[N]={};
	char id[N][10],name[N][100];
	double score[N][9];
	while (cin)
	{
		cin>>id[i]; 		
		cin >> name[i];
		
		for (j=0;j<8;j++)
		{
			cin >> score[i][j];
			s += score[i][j];
		}
		score[i][8]=s/8;
		
		s=0.0;
		i += 1;
	}
	
	for (t=0;t<i-1;t++)
	{
		r=0;
		while (flag[r] == 1)
			r+=1;
		
		max=score[r][8];

		for (p=0;p<i;p++)
		{
			if (score[p][8] >= max && flag[p] == 0)
			{
				max=score[p][8];
				imax=p;
			}
		}
		flag[imax]=1;
		cout << id[imax] << " ";
		cout << name[imax] << " ";
		cout << score[imax][8] << endl;
	}
	return 0;
}

