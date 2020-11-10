#include <iostream>
using namespace std;

bool zhishu(int i);

int main()
{
	int a=1,i;
	bool f;
	while ( a != 0)
	{
		f=true;
		a = 1;
		while (a % 2 !=0 || a <= 2)
		{
			cout << "请输入大于2偶数："; 
			cin >> a;		
		}
		for (i=2;i<=a/2;i++)
		{
			if (zhishu(i))
				if (zhishu(a-i))
					{
						cout << "对于输入的偶数，哥德巴赫猜想成立:" <<a<<"="<<i<<"+"<<a-i<<endl;
						f=false;
						break;
					}
		}
		if (f)
			cout <<"验证失败！"<<endl; 
	}
	return 0;
}

bool zhishu(int i)
{
	int k;
	bool flag=true;
	for (k=2;k<i/2;k++)
		if (i % k == 0)
		{
			flag = false;
			break;
		}	
	return flag;
}
