f_seq=open('EnzymeCut_SeqInfo.txt','r')
f_enzy=open('EnzymeCut_EnzyInfo.txt','r')
ls,le,c=[],[],0
for line in f_seq:
    line=line.replace('\n','')
    if c!=0:
        ls.append(line.upper())
    else:
        c=1
        print('Input Seq:')
        print(line)
s=''.join(ls)
print(s)

for line in f_enzy:
    line=line.replace('\n','')
    le.append(line.split(','))
print('\nInput Enzyme:')
for item in le:
    print(item[0],end=',')
    item[1]=item[1].upper()
print('\n')

site,ans=[],[]
for item in le:
    t,p=0,0
    while t!=-1:
        t=s.find(item[1],p)
        if t!=-1:
            p=t+eval(item[2])
            for i in site:
                if abs(p-i)<=12:    #设置两酶切位点距离阈值
                    break
            else:
                site.append(p)

p=0
site=sorted(site)
for i in site:
    ans.append([s[p:i],i-p])
    p=i
ans.append([s[p:],len(s)-p])

t=1
for item in ans:
    print('Frag',t,':',item[0],'Length:',item[1])
    t+=1

print('\nPress <Enter> to quit...',end='')
input()

