s=input('Input Seq:')
f=input('Input Frag:')
t=0
p=0
ans=[]
while t!=-1:
    t=s.find(f,p)
    p=t+1
    ans.append(str(p))
print('Position:',' '.join(ans[:-1]))

print('\nPress <Enter> to quit...',end='')
input()
