#include <iostream>
using namespace std;

struct node
{
	int content;
	node *next;
};
node *head=NULL;

node *mycreate(node *head);
node *mydelete(node *head);
void myprint(node *head);
void myfree(node *head);

int main()
{
	head=mycreate(head);
	//myprint(head);cout<<endl;
	head=mydelete(head);
	myprint(head);
	myfree(head);
}

node *mycreate(node *head)
{
	node *p=new node;
	head = p;
	p->next = NULL;
	cin >> p->content;
	node *q=p;
	while (p->content !=-1)
	{
		p=new node;
		cin >> p->content;
		q->next = p;
		p->next = NULL;
		q=p;
	}
	return head;
}

node *mydelete(node *head)
{
	node *h=head;
	while (h->content != -1)
	{
		node *p=h->next;
		while (p->content != -1)
			if (h->content == p->content)
			{ 
				h->next = p->next;
				node *q=p;
				p=p->next;
				delete q;
			} 
			else break;
		h = h->next;
	}
	return head;
}

void myprint(node *head)
{
	node *h=head;
	while (h->content !=-1)
	{
		cout << h->content <<" ";
		h = h->next;
	}
}

void myfree(node *head)
{
	while(head)
	{
		node *cur;
		cur = head;
		head = head->next;
		delete cur;
	}
}
