s=input()
d=eval(s)
s={}
try:
    l1=list(d.keys())
    l2=list(d.values())
    for i in range(len(l1)):
        s[l2[i]]=l1[i]
    print(s)
except:
    print('输入错误')
