#include <iostream>
using namespace std;

int main()
{
	int i=0,j,k;
	cin >> k;
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
	while (p->content != -1) //using -1 to end
	{
		p=new node;
		cin>>p->content;
		q->next=p;
		p->next=NULL;
		q=p;
		i++;
	}
	node *t=head;
	for (j=1;j<i-k;j++)
		t=t->next;
	cout << t->content;
	delete p;
	delete q;
	delete t;
	return 0;
} 
