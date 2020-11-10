#include<stdio.h>
int n,i,j,ans;

int main()
{
	while(scanf("%d",&n)!=EOF)
	{
		ans=0;
		for(i=5;i<=n;i*=5)
			for(j=1;j*i<=n;j++)
				ans++;
		printf("%d\n",ans);
	}
	return 0;
}
