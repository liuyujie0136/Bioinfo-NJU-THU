#include <iostream>
using namespace std;
#define N 100

int main()
{
	int i,j,k=0,t=0,m=0;
	int a1[N]={},a2[N]={},a3[N]={};
	bool flag=true;
	struct node
	{
		int content;
		node *next;
	};
	node *head=NULL;
	node *p=new node;
	head = p;
    cin>>p->content; 
	p->next=NULL;
	node *q=p;
	while (p->content != -2) //using -1 as the end of the first node &-2 as the end of the second node
	{
		if (p->content != -1 && flag)
			k++;
		if (p->content == -1)
			flag=false;
		p=new node;
		cin>>p->content;
		q->next=p;
		p->next=NULL;
		q=p;
		m++;
	}
	
	for (node *p=head; p!= NULL; p=p->next)
	{
		if (t<=k)
		{
			a1[t]=p->content;
			t++;
		}
		else
		{
			a2[t-k-1]=p->content;
			t++;
		}
	}
	
	for (i=0;i<k;i++)	//交集 
		for (j=0;j<m-k-1;j++)
			if (a1[i]==a2[j])
				cout << a1[i]<<" ";
	cout <<endl;
	
	for (i=0;i<k;i++)	//并集 
		a3[i]=a1[i];
	for (j=0;j<m-k-1;j++)
		a3[k+j]=a2[j];
	for (i=0;i<m-1;i++)
		for (j=i+1;j<m-1;j++)
			if (a3[i]==a3[j])
				a3[i]=-1;
	for (i=0;i<m-1;i++)
		if (a3[i] != -1)
			cout << a3[i] <<" ";
	cout <<endl;
	
	
	for (i=0;i<k;i++)	//差集 1-2 
	{
		flag=true;
		for (j=0;j<m-k-1;j++)
			if (a1[i]==a2[j])
				flag=false;
		if (flag)
			cout << a1[i]<<" ";
	}		
	cout <<endl;
	delete p;
	delete q;
	return 0;
}
