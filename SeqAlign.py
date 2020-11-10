def test(nt1,nt2):
    if nt1==nt2:
        return 2
    if nt1!=nt2:
        return -1
#计算得分

s1=input('Input Seq1:')
s2=input('Input Seq2:')
l1=len(s1)
l2=len(s2)
m=[]
for i in range((l1+1)*(l2+1)):
    m+='0'
k=0
for i in range(l1+1):
    m[i]=k
    k-=1
k=0
for i in range(0,(l1+1)*(l2+1),l1+1):
    m[i]=k
    k-=1
#初始化

for i in range(l2):
    for j in range(l1):
        m[(i+1)*(l1+1)+j+1]=max(m[i*(l1+1)+j]+test(s1[j],s2[i]),m[(i+1)*(l1+1)+j]-1,m[i*(l1+1)+j+1]-1)
#计算得分矩阵
'''
k=0
print('{:^12}'.format(' '),end=' ')
for i in range(l1):
    print('{:^6}'.format(s1[i]),end=' ')
print('\n\n','{:^5}'.format(' '),end=' ')
for j in range(l1+1):
    print('{:^5}'.format(m[k]),end=' ')
    k+=1
print('\n')
for i in range(l2):
    print('{:^5}'.format(s2[i]),end=' ')
    for j in range(l1+1):
        print('{:^5}'.format(m[k]),end=' ')
        k+=1
    print('\n')
#打印得分矩阵'''

ls1,ls2=list(s1),list(s2)
i,j=l1,l2
while i!=0 and j!=0:
    if s1[i-1]==s2[j-1] and i>0 and j>0:
            i-=1
            j-=1
    else:
        a=max(m[(j-1)*(l1+1)+i-1],m[j*(l1+1)+i-1],m[(j-1)*(l1+1)+i])
        if a==m[(j-1)*(l1+1)+i-1]:
            i-=1
            j-=1
        elif a==m[j*(l1+1)+i-1]:
            i-=1
            ls2.insert(j,'-')
        elif a==m[(j-1)*(l1+1)+i]:
            j-=1
            ls1.insert(i,'-')

if i==0 and s1[0]!=s2[0]:
    ls1.insert(0,'-')
if j==0 and s1[0]!=s2[0]:
    ls2.insert(0,'-')
#追溯序列

ss1,ss2=''.join(ls1),''.join(ls2)
print('\nAlignment:')
print('Seq1:',ss1)
print('Seq2:',ss2)
#打印序列

print('\nPress <Enter> to quit...',end='')
input()



