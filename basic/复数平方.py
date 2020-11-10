s=input()
f1=1
f2=0
f3=0
for item in s:
    if '0'<=item<='9':
        ans=''
    elif item=='+':
        f1=0
        f2=1
    elif item=='j':
        f1=0
        f3=1
    elif item=='.':
        f1=1
    else:
        f1=0
if f2==1 and f3==1:
    f1=1
if f1==0:
    ans='输入有误'
elif f2==1 and f3==1:
    t=s.split('+')
    a=eval(t[0])
    b=eval(t[1][:-1])
    aa=a*a-b*b
    bb=2*a*b
    ans='('+str(aa)+'+'+str(bb)+'j'+')'
else:
    ans=(eval(s))**2
print(ans)
